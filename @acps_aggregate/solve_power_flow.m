function [x, success, i] = solve_power_flow(obj, mpc, mpopt)
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

%% set up Newton solver options
opt = struct( ...
        'verbose',    mpopt.verbose, ...
        'tol',        mpopt.pf.tol, ...
        'max_it',     mpopt.pf.nr.max_it, ...
        'lin_solver', mpopt.pf.nr.lin_solver ...
    );

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(mpc.bus, mpc.gen);

%% create x0 for Newton power flow
va = obj.params_var('va');
vm = obj.params_var('vm');
zr = obj.params_var('zr');
zi = obj.params_var('zi');
v_ = vm .* exp(1j * va);
z_ = zr + 1j * zi;
x0 = [va([pv; pq]); vm(pq)];

fcn = @(x)power_flow_equations(obj, x, va, vm, z_, ref, pv, pq);

[x, success, i] = newton_solver(x0, fcn, opt);
