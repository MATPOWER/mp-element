classdef mpe_network_acps_test_nln < mpe_network_acp% & mp_model_acps

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
        function obj = mpe_network_acps_test_nln()
            obj@mpe_network_acp();
            obj.element_classes = ...
                { @mpe_bus_acp, @mpe_gen_acp_nln, @mpe_load_acp_nln, @mpe_branch_acp_nln, @mpe_shunt_acp_nln, @mpe_gizmo_acp_nln };

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  INIT_SET_TYPES() is called explicitly as needed
            %%              (if obj.node is empty) in CREATE_MODEL() and
            %%              DISPLAY(), after object construction, but before
            %%              object use.
        end


        %%-----  PF methods  -----
        function add_pf_vars(obj, nm, om, mpc, mpopt)
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
                    error('handling of indexed sets not implmented here (yet)');
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
                    error('handling of indexed sets not implmented here (yet)');
                end
            end
        end

        function [v_, z_] = pfx2vz(obj, x, ad)
            %% update v_, z_ from x
            ad.v1([ad.pv; ad.pq]) = x(1:ad.npv+ad.npq);         %% va
            ad.v2(ad.pq)          = x(ad.npv+ad.npq+1:end);     %% vm
            v_ = ad.v2 .* exp(1j * ad.v1);
            z_ = ad.zr + 1j * ad.zi;
        end

        function [f, J] = power_flow_equations(obj, x, ad)
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
            SS = C * S;
            f = [real(SS(pvq)); imag(SS(ad.pq))];
        end

        function add_pf_node_balance_constraints(obj, om, mpc, mpopt)
            %% power balance constraints
            ad = om.get_userdata('power_flow_aux_data');
            fcn = @(x)power_flow_equations(obj, x, ad);
            om.add_nln_constraint({'Pmis', 'Qmis'}, [ad.npv+ad.npq;ad.npq], 1, fcn, []);
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
