classdef mp_network_accs < mp_network_acc% & mp_form_accs

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
        function pf_add_vars(obj, mm, nm, dm, mpopt)
            %% get model variables
            vvars = obj.model_vvars();

            %% index vectors
            ad = mm.get_userdata('aux_data');
            pqv = [ad.pq; ad.pv];

            %% voltage real part
            st = obj.(vvars{1});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    mm.add_var(name, ad.npq+ad.npv, d.v0.(name)(pqv), d.vl.(name)(pqv), d.vu.(name)(pqv));
                else
                    error('mp_network_accs/pf_add_vars: handling of indexed sets not implmented here (yet)');
                end
            end

            %% voltage imaginary part
            st = obj.(vvars{2});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    mm.add_var(name, ad.npq+ad.npv, d.v0.(name)(pqv), d.vl.(name)(pqv), d.vu.(name)(pqv));
                else
                    error('mp_network_accs/pf_add_vars: handling of indexed sets not implmented here (yet)');
                end
            end
        end

        function [vx_, z_] = pf_convert_x(obj, x, ad)
            %% update v_, z_ from x
            pqv = [ad.pq; ad.pv];
            ad.v1(pqv) = x(1:ad.npv+ad.npq);                    %% vr
            ad.v2(pqv) = x(ad.npv+ad.npq+1:2*ad.npv+2*ad.npq);  %% vi
            vx_ = ad.v1 + 1j * ad.v2;
            z_ = ad.zr + 1j * ad.zi;
            if nargout < 2
                vx_ = [vx_; z_];
            end
        end

        function [f, J] = pf_node_balance_equations(obj, x, ad)
            %% index vector
            pqv = [ad.pq; ad.pv];

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.pf_convert_x(x, ad);

            %% incidence matrix
            C = obj.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                [S, Svr, Svi] = obj.port_inj_power([v_; z_], 1);
                dV2_dVr = sparse(1:ad.npv, ad.npq+(1:ad.npv), 2*real(v_(ad.pv)), ad.npv, ad.npv+ad.npq);
                dV2_dVi = sparse(1:ad.npv, ad.npq+(1:ad.npv), 2*imag(v_(ad.pv)), ad.npv, ad.npv+ad.npq);

                SSvr = C * Svr;
                SSvi = C * Svi;
                J = [   real(SSvr(pqv,  pqv))   real(SSvi(pqv,  pqv));
                        imag(SSvr(ad.pq, pqv))  imag(SSvi(ad.pq, pqv));
                        dV2_dVr dV2_dVi ];
            else
                %% get port power injections (w/o derivatives)
                S = obj.port_inj_power([v_; z_], 1);
            end

            %% nodal power balance + voltage magnitude mismatch
            SS = C * S;
            vmm = v_(ad.pv) .* conj(v_(ad.pv)) - ad.v1(ad.pv).^2 - ad.v2(ad.pv).^2;
            f = [real(SS(pqv)); imag(SS(ad.pq)); vmm];
        end

        function pf_add_node_balance_constraints(obj, mm, dm, mpopt)
            %% power balance constraints
            ad = mm.get_userdata('aux_data');
            fcn = @(x)pf_node_balance_equations(obj, x, ad);
            mm.add_nln_constraint({'Pmis', 'Qmis', 'Vmis'}, [ad.npv+ad.npq;ad.npq;ad.npv], 1, fcn, []);
        end


        %%-----  OPF methods  -----
        function opf_add_node_balance_constraints(obj, mm)
            %% power balance constraints
            nn = obj.node.N;            %% number of nodes
            fcn_mis = @(x)opf_power_balance_fcn(obj, obj.opf_convert_x(x));
            hess_mis = @(x, lam)opf_power_balance_hess(obj, ...
                obj.opf_convert_x(x), lam);
            mm.add_nln_constraint({'Pmis', 'Qmis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end
    end     %% methods
end         %% classdef
