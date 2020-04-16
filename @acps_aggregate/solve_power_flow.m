function [v_, success, i, data] = solve_power_flow(obj, mpc, mpopt)
%SOLVE_POWER_FLOW  Solves Newton power flow
%   SUCCESS = SOLVE_POWER_FLOW(OBJ, MPC)
%
%   Inputs
%       OBJ - 
%       MPC - 
%       MPOPT - 
%
%   Returns
%       SUCCESS - 
%
%   See also ...

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% MATPOWER options
if nargin < 3
    mpopt = mpoption;
end

if mpopt.verbose, fprintf('-----  solve_power_flow()  -----\n'); end

%% set up Newton solver options
opt = struct( ...
        'verbose',    mpopt.verbose, ...
        'tol',        mpopt.pf.tol, ...
        'max_it',     mpopt.pf.nr.max_it, ...
        'lin_solver', mpopt.pf.nr.lin_solver ...
    );

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(mpc.bus, mpc.gen);
node_types = struct('ref', ref, 'pv', pv, 'pq', pq, ...
        'nref', length(ref), 'npv', length(pv), 'npq', length(pq));

% keyboard

%% create x0 for Newton power flow
vvars = obj.model_vvars();
zvars = obj.model_zvars();
v1 = obj.params_var(vvars{1});
v2 = obj.params_var(vvars{2});
zr = obj.params_var(zvars{1});
zi = obj.params_var(zvars{2});

x0 = obj.vz2pfx(v1, v2, zr, zi, node_types);

fcn = @(x)power_flow_equations(obj, x, v1, v2, zr, zi, node_types);

[x, success, i] = newton_solver(x0, fcn, opt);

[v_, z_] = pfx2vz(obj, x, v1, v2, zr, zi, node_types);
