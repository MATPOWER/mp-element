classdef mp_math_pf_acps < mp_math_pf & mm_pf_shared_acps
%MP_MATH_PF_ACPS  MATPOWER mathematical model for AC power flow (PF) problem.
%   ?
%
%   MP_MATH_PF_ACPS ... power flow ...
%
%   Properties
%       ? - ?
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        %% constructor
        function obj = mp_math_pf_acps()
            obj@mp_math_pf();
            obj.element_classes = { @mme_buslink_pf_acp };
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            alg = mpopt.pf.alg;
            ad = obj.aux_data;
            
            %% power balance constraints
            switch alg
                case  {'FDXB', 'FDBX'}
                    fcn = @(x)pf_node_balance_equations(obj, x, nm, ad, 1);
                otherwise
                    fcn = @(x)pf_node_balance_equations(obj, x, nm, ad);
            end
            obj.add_nln_constraint({'Pmis', 'Qmis'}, [ad.npv+ad.npq;ad.npq], 1, fcn, []);
        end

        function x = gs_x_update(obj, x, f, nm, dm, mpopt);
            alg = mpopt.pf.alg;
            ad = obj.aux_data;

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.pf_convert_x(x, nm, ad, 1);

            [pv, pq, npv, npq, Y] = deal(ad.pv, ad.pq, ad.npv, ad.npq, ad.Y);
            
            %% total nodal complex bus power extractions
            SS = zeros(size(v_));
            SS(pv) = f(1:npv);
            SS(pq) = f(npv+1:npv+npq) + 1j * f(npv+npq+1:npv+2*npq);
%            SS = C * nm.port_inj_power([v_; z_], 1);

            %% complex net nodal injection (from all but constant Z elements)
            S0 = v_ .* conj(Y * v_) - SS;

            %% update voltage
            %% at PQ buses
            for k = pq'
                v_(k) = v_(k) + (conj(S0(k)/v_(k)) - Y(k,:) * v_) / Y(k, k);
            end

            %% at PV buses
            if npv
                for k = pv'
                    S0(k) = real(S0(k)) + 1j * imag( v_(k) * conj(Y(k,:) * v_) );
                    v_(k) = v_(k) + (conj(S0(k)/v_(k)) - Y(k,:) * v_) / Y(k, k);
                end
            end

            x = [angle(v_([pv; pq])); abs(v_(pq))];
        end

        function x = zg_x_update(obj, x, f, nm, dm, mpopt);
            alg = mpopt.pf.alg;
            ad = obj.aux_data;

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.pf_convert_x(x, nm, ad, 1);

            [pv, pq, ref, npv, npq] = deal(ad.pv, ad.pq, ad.ref, ad.npv, ad.npq);
            pvq = [pv; pq];

            %% build and cache S0 and Y21_v1
            if isfield(ad, 'S0')
                [Y21_v1, S0] = deal(ad.Y21_v1, ad.S0);
            else
                Y21_v1 = ad.Y21 * v_(ref);
                SS = zeros(size(v_));
                SS(pv) = f(1:npv);
                SS(pq) = f(npv+1:npv+npq) + 1j * f(npv+npq+1:npv+2*npq);
    %            SS = C * nm.port_inj_power([v_; z_], 1);

                %% complex net nodal injection (from all but constant Z elements)
                S0 = v_ .* conj(ad.Y * v_) - SS;
                S0(ref) = 0;

                %% cache 'em
                [ad.Y21_v1, ad.S0] = deal(Y21_v1, S0);
                obj.aux_data = ad;
            end

            if npv  %% update Q injections at PV buses based on vm mismatch
                %% compute and cache initial vm at PV buses and
                %% factored fast-decoupled Bpp matrix
                if isfield(ad, 'vmpv0')
                    [vmpv, vmpv0, Bpp, LBpp, UBpp, pBpp, iqBpp] = ...
                        deal(ad.vmpv, ad.vmpv0, ad.Bpp, ad.LBpp, ad.UBpp, ad.pBpp, ad.iqBpp);
                else
                    vmpv0 = abs(v_(pv));
                    vmpv = vmpv0;

                    %% modify data model to form Bpp (B double prime)
                    dm2 = dm.fdpf_B_matrix_models('FDBX');

                    %% build network models and get admittance matrices
                    nm2 = feval(class(nm)).build(dm2);
                    [Y2, L, M] = nm2.get_params([], {'Y', 'L', 'M'});
                    if any(any(L)) || any(any(M))
                        error('mp_math_pf_acps/zg_x_update: B matrix for Z-bus Gauss w/PV buses not implemented for models with non-zero L and/or M matrices.')
                    end
                    Bpp = -nm2.C * imag(Y2) * nm2.C';

                    [LBpp, UBpp, pBpp, qBpp] = lu(Bpp(pq, pq), 'vector');
                    [junk, iqBpp] = sort(qBpp);

                    %% cache 'em
                    [ad.vmpv, ad.vmpv0, ad.Bpp, ad.LBpp, ad.UBpp, ad.pBpp, ad.iqBpp] = ...
                        deal(vmpv, vmpv0, Bpp, LBpp, UBpp, pBpp, iqBpp);
                    obj.aux_data = ad;
                end

                %% compute voltage mismatches at PV buses
                v_(pv) = vmpv .* v_(pv) ./ abs(v_(pv));
                dV = vmpv0 - vmpv;
%                 [max_dV, k] = max(abs(dV));
%                 fprintf('       %10.3e', max_dV)

                %% compute Q injection at current V
                %% (sometimes improves convergence)
                Qpv = imag( v_(pv) .* conj(ad.Y(pv, :) * v_) );
                S0(pv) = S0(pv) + 1j * (Qpv - imag(S0(pv)));

                % dVpq = Bpp(pq, pq) \ (-Bpp(pq, pv) * dV);
                dVpq = UBpp \  (LBpp \ (-Bpp(pq(pBpp), pv) * dV));
                dVpq = dVpq(iqBpp);
                dQ = Bpp(pv, pq) * dVpq + Bpp(pv, pv) * dV;

                %% update S0
                S0(pv) = S0(pv) + 1j * dQ;
            end

            %% complex current injections
            I2 = conj(S0(pvq) ./ v_(pvq));

            V2 = ad.U \  (ad.L \ (I2(ad.p) - Y21_v1(ad.p)));
            V2 = V2(ad.iq);

            v_(pv) = V2(1:npv);
            v_(pq) = V2(npv+1:npv+npq);
            obj.aux_data.vmpv = abs(v_(pv));

            x(1:ad.npv+ad.npq) = angle(v_(pvq));
            x(ad.npv+ad.npq+1:ad.npv+2*ad.npq) = abs(v_(pq));
        end
    end     %% methods
end         %% classdef
