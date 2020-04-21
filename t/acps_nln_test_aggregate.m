classdef acps_nln_test_aggregate < acp_aggregate% & acps_model

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
        function obj = acps_nln_test_aggregate()
            obj@acp_aggregate();
            obj.element_classes = ...
                { @acp_bus, @acp_nln_gen, @acp_nln_load, @acp_nln_branch, @acp_nln_shunt, @acp_nln_gizmo };
            if isempty(obj.node)    %% skip if constructed from existing object
                obj.init_set_types();   %% should be called in mp_idx_manager
                                        %% constructor, if not for:
            end                         %% https://savannah.gnu.org/bugs/?52614
        end


        %%-----  PF methods  -----
        function x = vz2pfx(obj, va, vm, zr, zi, t, ad)
            %% update x from va, vm, zr, zi
            x = [va([t.pv; t.pq]); vm(t.pq)];
        end

        function [v_, z_] = pfx2vz(obj, x, va, vm, zr, zi, t, ad)
            %% update v_, z_ from x
            va([t.pv; t.pq]) = x(1:t.npv+t.npq);
            vm(t.pq) = x(t.npv+t.npq+1:end);
            v_ = vm .* exp(1j * va);
            z_ = zr + 1j * zi;
        end

        function [F, J] = power_flow_equations(obj, x, va, vm, zr, zi, t, ad)
            %% index vector
            pvq = [t.pv; t.pq];

            %% update model state ([v_; z_]) from power flow state (x)
            [v_, z_] = pfx2vz(obj, x, va, vm, zr, zi, t);

            %% incidence matrix
            C = obj.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                [S, Sva, Svm] = obj.port_inj_power([v_; z_], 1);

                SSva = C * Sva;
                SSvm = C * Svm;
                J = [   real(SSva(pvq,  pvq)) real(SSvm(pvq,  t.pq));
                        imag(SSva(t.pq, pvq)) imag(SSvm(t.pq, t.pq))  ];
            else
                %% get port power injections (w/o derivatives)
                S = obj.port_inj_power([v_; z_], 1);
            end

            %% nodal power balance
            SS = C * S;
            F = [real(SS(pvq)); imag(SS(t.pq))];
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
