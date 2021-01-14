classdef mp_network_acpi < mp_network_acp & mp_form_acpi

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
        function ad = pf_aux_data(obj, dm, mpopt)
            %% call parent method
            ad = pf_aux_data@mp_network_ac(obj, dm, mpopt);

            %% build additional aux data
            N = obj.C(ad.pv, :) * obj.N;%% z coefficients for z @ PV nodes
            [ii, jj, ss] = find(N);     %% deconstruct and recreate with
            [~, ia] = unique(ii, 'first');  %% only 1st non-zero in each row
            N = sparse(ii(ia), jj(ia), ss(ia), ad.npv, size(N, 2)) * obj.D';
            k = find(any(N, 1));        %% indices of PV node z-vars

            %% save additional aux data
            ad.N = N(:,k);          %% coefficients for zi(k) corrsp to PV nodes
            ad.invN = inv(ad.N);    %% inverse of N, typically equal to N = -eye()
            ad.k = k;               %% indices of PV node z-vars
        end

        function obj = pf_add_vars(obj, mm, nm, dm, mpopt)
            %% get model variables
            vvars = obj.model_vvars();

            %% index vectors
            ad = mm.get_userdata('aux_data');
            pvq = [ad.pv; ad.pq];

            %% voltage angles
            st = obj.(vvars{1});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    mm.add_var(name, ad.npv+ad.npq, d.v0.(name)(pvq), d.vl.(name)(pvq), d.vu.(name)(pvq));
                else
                    error('mp_network_acpi/pf_add_vars: handling of indexed sets not implmented here (yet)');
                end
            end

            %% reactive injections
            v_ = ad.v1 + 1j * ad.v2;
            z_ = ad.zr + 1j * ad.zi;
            Qpv = obj.C(ad.pv, :) * imag( obj.port_inj_power([v_; z_], 1) );
            Qg_pv = Qpv - ad.N * ad.zi(ad.k);
            mm.add_var('Qg_pv', ad.npv, Qg_pv);

            %% voltage magnitudes
            st = obj.(vvars{2});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    mm.add_var(name, ad.npq, d.v0.(name)(ad.pq), d.vl.(name)(ad.pq), d.vu.(name)(ad.pq));
                else
                    error('mp_network_acpi/pf_add_vars: handling of indexed sets not implmented here (yet)');
                end
            end
        end

        function [vx_, z_, x_] = pf_convert_x(obj, mmx, ad, only_v)
            %% x = obj.pf_convert(mmx, ad)
            %% [v, z] = obj.pf_convert(mmx, ad)
            %% [v, z, x] = obj.pf_convert(mmx, ad)
            %% ... = obj.pf_convert(mmx, ad, only_v)

            %% update v_, z_ from mmx
            iN = ad.npv + ad.npq;           ad.v1([ad.pv; ad.pq]) = mmx(1:iN);  %% va
            i1 = iN+1;  iN = iN + ad.npv;   Qg_pv = mmx(i1:iN);                 %% Qg_pv
            i1 = iN+1;  iN = iN + ad.npq;   ad.v2(ad.pq) = mmx(i1:iN);          %% vm
            vx_ = ad.v2 .* exp(1j * ad.v1);
            ad.zi(ad.k) = -ad.N \ Qg_pv;
            z_ = ad.zr + 1j * ad.zi;

            %% update z, if requested
            if nargin < 4 || ~only_v
                z_ = obj.pf_update_z(vx_, z_, ad);
            end

            if nargout < 2
                vx_ = [vx_; z_];
            elseif nargout > 2
                x_ = [vx_; z_];
            end
        end

        function [f, J] = pf_node_balance_equations(obj, x, ad)
            %% index vector
            pvq = [ad.pv; ad.pq];

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.pf_convert_x(x, ad, 1);

            %% incidence matrix
            C = obj.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                [I, Iva, Ivm, Izr, Izi] = obj.port_inj_current([v_; z_], 1);

                IIva = C * Iva;
                IIvm = C * Ivm;
                IIzi = C * Izi;
                IIvm(:, ad.pv) = -IIzi(:, ad.k) * ad.invN;  %% dImis_dQg

                J = [   real(IIva(pvq, pvq))    real(IIvm(pvq, pvq));
                        imag(IIva(pvq, pvq))    imag(IIvm(pvq, pvq))    ];
            else
                %% get port power injections (w/o derivatives)
                I = obj.port_inj_current([v_; z_], 1);
            end

            %% nodal power balance
            II = C * I;
            f = [real(II(pvq)); imag(II(pvq))];
        end

        function obj = pf_add_node_balance_constraints(obj, mm, dm, mpopt)
            %% power balance constraints
            ad = mm.get_userdata('aux_data');
            npvq = ad.npv+ad.npq;
            fcn = @(x)pf_node_balance_equations(obj, x, ad);
            mm.add_nln_constraint({'Irmis', 'Iimis'}, [npvq;npvq], 1, fcn, []);
        end


        %%-----  OPF methods  -----
        function opf_add_node_balance_constraints(obj, mm)
            %% power balance constraints
            nn = obj.node.N;            %% number of nodes
            fcn_mis = @(x)opf_current_balance_fcn(obj, obj.opf_convert_x(x));
            hess_mis = @(x, lam)opf_current_balance_hess(obj, ...
                obj.opf_convert_x(x), lam);
            mm.add_nln_constraint({'rImis', 'iImis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end

        function [lamP, lamQ] = opf_node_power_balance_prices(obj, mm)
            %% shadow prices on node power balance
            nne = mm.get_idx('nle');

            %% convert current balance shadow prices to equivalent lamP and lamQ
            %% P + jQ = (Vr + jVi) * (M - jN)
            %% M = (Vr P + Vi Q) / (Vr^2 + Vi^2)
            %% N = (Vi P - Vr Q) / (Vr^2 + Vi^2)
            %% lamP = df/dP = df/dM * dM/dP + df/dN + dN/dP
            %% lamQ = df/dQ = df/dM * dM/dQ + df/dN + dN/dQ
            V = obj.soln.v;
            lambda = mm.soln.lambda;
            VV = V ./ (V .* conj(V));   %% V / Vm^2
            VVr = real(VV);
            VVi = imag(VV);
            lamM = lambda.eqnonlin(nne.i1.rImis:nne.iN.rImis);
            lamN = lambda.eqnonlin(nne.i1.iImis:nne.iN.iImis);
            lamP = (VVr.*lamM + VVi.*lamN);
            lamQ = (VVi.*lamM - VVr.*lamN);
        end
    end     %% methods
end         %% classdef
