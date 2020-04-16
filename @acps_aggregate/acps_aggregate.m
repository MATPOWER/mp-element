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
        function add_opf_node_balance_constraints(obj, om)
            %% power balance constraints
            nn = obj.node.N;            %% number of nodes
            fcn_mis = @(x)opf_power_balance_fcn(obj, x);
            hess_mis = @(x, lam)opf_power_balance_hess(obj, x, lam);
            om.add_nln_constraint({'Pmis', 'Qmis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end

        [F, J] = power_flow_equations(obj, x, va, vm, z_, ref, pv, pq)

        [x, success, i] = solve_power_flow(obj, mpc, mpopt)

        function x = vz2pfx(obj, va, vm, zr, zi, t)
            x = [va([t.pv; t.pq]); vm(t.pq)];
        end

        function [v_, z_] = pfx2vz(obj, x, va, vm, zr, zi, t)
            va([t.pv; t.pq]) = x(1:t.npv+t.npq);
            vm(t.pq) = x(t.npv+t.npq+1:end);
            v_ = vm .* exp(1j * va);
            z_ = zr + 1j * zi;
        end
    end     %% methods
end         %% classdef
