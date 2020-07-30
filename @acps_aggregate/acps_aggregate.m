classdef acps_aggregate < acp_aggregate% & acps_model

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
        function add_pf_vars(obj, asm, om, ad, mpc, mpopt)
            %% get model variables
            vvars = obj.model_vvars();

            %% index vectors
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
            SS = C * S;
            f = [real(SS(pvq)); imag(SS(ad.pq))];
        end

        function add_pf_node_balance_constraints(obj, om, ad)
            %% power balance constraints
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
