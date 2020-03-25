function [x, success, i] = solve_opf(obj, mpc, mpopt)
%SOLVE_OPF  Solves AC OPF
%   SUCCESS = SOLVE_OPF(OBJ, MPC)
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

%% create opf_model object
om = obj.setup_opf(mpc, mpopt);

%% solve it
opt = mpopt2nlpopt(mpopt, om.problem_type(), 'DEFAULT');
[x, f, eflag, output, lambda] = om.solve(opt);
success = (eflag > 0);

% om
% vv = om.get_idx('var');
% x(vv.i1.Pg:vv.iN.Pg)
% keyboard

i = output.iterations;
