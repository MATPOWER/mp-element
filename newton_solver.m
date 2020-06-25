function [x, success, i] = newton_solver(x0, fcn, opt)
%NEWTON_SOLVER  Solves F(x) = 0 via Newton's method
%   [X, SUCCESS, I] = NEWTON_SOLVER(X0, FCN)
%   [X, SUCCESS, I] = NEWTON_SOLVER(X0, FCN, OPT)
%
%       X0 - full system admittance matrix (for all buses)
%       FCN - handle to function that computes F(x) and it's Jacobian
%           F = FCN(X)
%           [F, J] = FCN(X)
%       OPT - (optional) option struct, with fields
%           verbose - 
%           max_it - 
%           tol - 
%           lin_solver - 
%
%   Returns
%       X - solution
%       SUCCESS - 1 if converged successfully, 0 otherwise
%       I - number of iterations performed
%
%   See also ...

%   MATPOWER
%   Copyright (c) 1996-2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% default arguments
if nargin < 3
    opt = struct( ...
            'verbose', 2, ...
            'max_it', 10, ...
            'tol', 1e-8, ...
            'lin_solver', '\' ...
        );
end

%% options
tol         = opt.tol;
max_it      = opt.max_it;
lin_solver  = opt.lin_solver;

%% initialize
success = 0;
i = 0;
x = x0;

%% evaluate F(x0)
[F, J] = fcn(x);

%% check tolerance
normF = norm(F, inf);
if opt.verbose > 1
    fprintf('\n it     max residual');
    fprintf('\n----  ----------------');
    fprintf('\n%3d     %10.3e', i, normF);
end
if normF < tol
    success = 1;
    if opt.verbose > 1
        fprintf('\nConverged!\n');
    end
end

% %% attempt to pick fastest linear solver, if not specified
% if isempty(lin_solver)
%     nx = length(F);
%     if nx <= 10 || have_fcn('octave')
%         lin_solver = '\';       %% default \ operator
%     else    %% MATLAB and nx > 10
%         lin_solver = 'LU3';     %% LU decomp with 3 output args, AMD ordering
%     end
% end

%% do Newton iterations
while (~success && i < max_it)
    %% update iteration counter
    i = i + 1;

    %% compute update step
    dx = J \ -F;
%     dx = mplinsolve(J, -F, lin_solver);

    %% update voltage
    x = x + dx;

    %% evalute F(x) and J(x)
    [F, J] = fcn(x);

    %% check for convergence
    normF = norm(F, inf);
    if opt.verbose > 1
        fprintf('\n%3d     %10.3e', i, normF);
    end
    if normF < tol
        success = 1;
        if opt.verbose
            fprintf('\nNewton''s method converged in %d iterations.\n', i);
        end
    end
end

if opt.verbose
    if ~success
        fprintf('\nNewton''s method did not converge in %d iterations.\n', i);
    end
end
