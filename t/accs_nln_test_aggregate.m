classdef accs_nln_test_aggregate < acc_aggregate% & accs_model

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
        function obj = accs_nln_test_aggregate()
            obj@acc_aggregate();
            obj.element_classes = ...
                { @acc_bus, @acc_nln_gen, @acc_nln_load, @acc_nln_branch, @acc_nln_shunt, @acc_nln_gizmo };
            if isempty(obj.node)    %% skip if constructed from existing object
                obj.init_set_types();   %% should be called in mp_idx_manager
                                        %% constructor, if not for:
            end                         %% https://savannah.gnu.org/bugs/?52614
        end


        %%-----  PF methods  -----
        function add_pf_vars(obj, asm, om, ad, mpc, mpopt)
            %% get model variables
            vvars = obj.model_vvars();

            %% index vectors
            pqv = [ad.pq; ad.pv];

            %% voltage real part
            st = obj.(vvars{1});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    om.add_var(name, ad.npq+ad.npv, d.v0.(name)(pqv), d.vl.(name)(pqv), d.vu.(name)(pqv));
                else
                    error('handling of indexed sets not implmented here (yet)');
                end
            end

            %% voltage imaginary part
            st = obj.(vvars{2});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    om.add_var(name, ad.npq+ad.npv, d.v0.(name)(pqv), d.vl.(name)(pqv), d.vu.(name)(pqv));
                else
                    error('handling of indexed sets not implmented here (yet)');
                end
            end
        end

        function [v_, z_] = pfx2vz(obj, x, ad)
            %% update v_, z_ from x
            pqv = [ad.pq; ad.pv];
            ad.v1(pqv) = x(1:ad.npv+ad.npq);        %% vr
            ad.v2(pqv) = x(ad.npv+ad.npq+1:end);    %% vi
            v_ = ad.v1 + 1j * ad.v2;
            z_ = ad.zr + 1j * ad.zi;
        end

        function [F, J] = power_flow_equations(obj, x, ad)
            %% index vector
            pqv = [ad.pq; ad.pv];

            %% update model state ([v_; z_]) from power flow state (x)
            [v_, z_] = obj.pfx2vz(x, ad);

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
            F = [real(SS(pqv)); imag(SS(ad.pq)); vmm];
        end

        function add_pf_node_balance_constraints(obj, om, ad)
            %% power balance constraints
            fcn = @(x)power_flow_equations(obj, x, ad);
            om.add_nln_constraint({'Pmis', 'Qmis', 'Vmis'}, [ad.npv+ad.npq;ad.npq;ad.npv], 1, fcn, []);
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
