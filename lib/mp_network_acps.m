classdef mp_network_acps < mp_network_acp% & mp_form_acps

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         va = [];
%         vm = [];
%     end
    
    methods
        %%-----  PF methods  -----
        function ad = power_flow_aux_data(obj, dm, mpopt)
            %% call parent method
            ad = power_flow_aux_data@mp_network_ac(obj, dm, mpopt);

            switch mpopt.pf.alg
                case 'GS'
                    ad.Y = obj.C * obj.get_params([], 'Y') * obj.C';
                case 'ZG'
                    pvq = [ad.pv; ad.pq];
                    Y = obj.C * obj.get_params([], 'Y') * obj.C';
                    Y21 = Y(pvq, ad.ref);
                    Y22 = Y(pvq, pvq);
                    [L, U, p, q] = lu(Y22, 'vector');
                    [junk, iq] = sort(q);
                    [ad.Y, ad.Y21, ad.L, ad.U, ad.p, ad.iq] = ...
                        deal(Y, Y21, L, U, p, iq);
            end
        end

        function add_pf_vars(obj, nm, om, dm, mpopt)
            %% get model variables
            vvars = obj.model_vvars();

            %% index vectors
            ad = om.get_userdata('power_flow_aux_data');
            pvq = [ad.pv; ad.pq];

            %% voltage angles
            st = obj.(vvars{1});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    om.add_var(name, ad.npv+ad.npq, d.v0.(name)(pvq), d.vl.(name)(pvq), d.vu.(name)(pvq));
                else
                    error('mp_network_acps/add_pf_vars: handling of indexed sets not implmented here (yet)');
                end
            end

            %% voltage magnitudes
            st = obj.(vvars{2});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    om.add_var(name, ad.npq, d.v0.(name)(ad.pq), d.vl.(name)(ad.pq), d.vu.(name)(ad.pq));
                else
                    error('mp_network_acps/add_pf_vars: handling of indexed sets not implmented here (yet)');
                end
            end
        end

        function [v_, z_] = pfx2vz(obj, x, ad)
            %% update v_, z_ from x
            ad.v1([ad.pv; ad.pq]) = x(1:ad.npv+ad.npq);                 %% va
            ad.v2(ad.pq)          = x(ad.npv+ad.npq+1:ad.npv+2*ad.npq); %% vm
            v_ = ad.v2 .* exp(1j * ad.v1);
            z_ = ad.zr + 1j * ad.zi;
        end

        function [f, J] = power_flow_equations(obj, x, ad, fdpf)
            %% index vector
            pvq = [ad.pv; ad.pq];

            %% update model state ([v_; z_]) from power flow state (x)
            [v_, z_] = obj.pfx2vz(x, ad);

            %% incidence matrix
            C = obj.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                [S, Sva, Svm] = obj.port_inj_power([v_; z_], 1);

                SSva = C * Sva;
                SSvm = C * Svm;
                J = [   real(SSva(pvq,   pvq))  real(SSvm(pvq,   ad.pq));
                        imag(SSva(ad.pq, pvq))  imag(SSvm(ad.pq, ad.pq))    ];
            else
                %% get port power injections (w/o derivatives)
                S = obj.port_inj_power([v_; z_], 1);
            end

            %% nodal power balance
            if nargin > 3 && fdpf
                SS = C * S ./ abs(v_);  %% for fast-decoupled formulation
            else
                SS = C * S;
            end
            f = [real(SS(pvq)); imag(SS(ad.pq))];
        end

        function add_pf_node_balance_constraints(obj, om, dm, mpopt)
            alg = mpopt.pf.alg;
            ad = om.get_userdata('power_flow_aux_data');
            
            %% power balance constraints
            switch alg
                case  {'FDXB', 'FDBX'}
                    fcn = @(x)power_flow_equations(obj, x, ad, 1);
                otherwise
                    fcn = @(x)power_flow_equations(obj, x, ad);
            end
            om.add_nln_constraint({'Pmis', 'Qmis'}, [ad.npv+ad.npq;ad.npq], 1, fcn, []);
        end

        function JJ = fd_jac_approx(obj, om, dm, mpopt)
            alg = mpopt.pf.alg;

            %% create copies of data model for building B prime, B double prime
            [dm1, dm2] = dm.fdpf_B_matrix_models(alg);

            %% build network models and get admittance matrices
            nm1 = feval(class(obj)).build(dm1, mpopt);
            nm2 = feval(class(obj)).build(dm2, mpopt);
            [Y1, L, M] = nm1.get_params([], {'Y', 'L', 'M'});
            Y2 = nm2.get_params();
            if any(any(L)) || any(any(M))
                error('mp_network_acps/df_jac_approx: fast-decoupled Jacobian approximation not implemented for models with non-zero L and/or M matrices.')
            end

            %% form reduced Bp and Bpp matrices
            ad = om.get_userdata('power_flow_aux_data');
            Cp  = nm1.C([ad.pv; ad.pq], :);
            Cpp = nm2.C(ad.pq, :);
            Bp  = -imag( Cp  * Y1 * Cp' );
            Bpp = -imag( Cpp * Y2 * Cpp' );
            JJ = {Bp, Bpp};
        end

        function x = gs_x_update(obj, x, f, om, dm, mpopt);
            alg = mpopt.pf.alg;
            ad = om.get_userdata('power_flow_aux_data');

            %% get model state ([v_; z_]) from power flow state (x)
            [v_, z_] = obj.pfx2vz(x, ad);

            [pv, pq, npv, npq, Y] = deal(ad.pv, ad.pq, ad.npv, ad.npq, ad.Y);
            
            %% total nodal complex bus power extractions
            SS = zeros(size(v_));
            SS(pv) = f(1:npv);
            SS(pq) = f(npv+1:npv+npq) + 1j * f(npv+npq+1:npv+2*npq);
%            SS = C * obj.port_inj_power([v_; z_], 1);

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


        function x = zg_x_update(obj, x, f, om, dm, mpopt);
            alg = mpopt.pf.alg;
            ad = om.get_userdata('power_flow_aux_data');

            %% get model state ([v_; z_]) from power flow state (x)
            [v_, z_] = obj.pfx2vz(x, ad);

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
    %            SS = C * obj.port_inj_power([v_; z_], 1);

                %% complex net nodal injection (from all but constant Z elements)
                S0 = v_ .* conj(ad.Y * v_) - SS;
                S0(ref) = 0;

                %% cache 'em
                [ad.Y21_v1, ad.S0] = deal(Y21_v1, S0);
                om.userdata.power_flow_aux_data = ad;
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
                    nm = feval(class(obj)).build(dm2, mpopt);
                    [Y2, L, M] = nm.get_params([], {'Y', 'L', 'M'});
                    if any(any(L)) || any(any(M))
                        error('mp_network_acps/zg_x_update: B matrix for Z-bus Gauss w/PV buses not implemented for models with non-zero L and/or M matrices.')
                    end
                    Bpp = -nm.C * imag(Y2) * nm.C';

                    [LBpp, UBpp, pBpp, qBpp] = lu(Bpp(pq, pq), 'vector');
                    [junk, iqBpp] = sort(qBpp);

                    %% cache 'em
                    [ad.vmpv, ad.vmpv0, ad.Bpp, ad.LBpp, ad.UBpp, ad.pBpp, ad.iqBpp] = ...
                        deal(vmpv, vmpv0, Bpp, LBpp, UBpp, pBpp, iqBpp);
                    om.userdata.power_flow_aux_data = ad;
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
                ad.S0 = S0;
            end

            %% complex current injections
            I2 = conj(S0(pvq) ./ v_(pvq));

            V2 = ad.U \  (ad.L \ (I2(ad.p) - Y21_v1(ad.p)));
            V2 = V2(ad.iq);

            v_(pv) = V2(1:npv);
            v_(pq) = V2(npv+1:npv+npq);
            om.userdata.power_flow_aux_data.vmpv = abs(v_(pv));

            x = [angle(v_(pvq)); abs(v_(pq))];
        end


        %%-----  OPF methods  -----
        function add_opf_node_balance_constraints(obj, om)
            %% power balance constraints
            nn = obj.node.N;            %% number of nodes
            fcn_mis = @(x)opf_power_balance_fcn(obj, x);
            hess_mis = @(x, lam)opf_power_balance_hess(obj, x, lam);
            om.add_nln_constraint({'Pmis', 'Qmis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end
    end     %% methods
end         %% classdef
