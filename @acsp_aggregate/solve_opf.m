function [x, success, i] = solve_opf(obj, mpc, mpopt)
%SOLVE_OPF  Solves Newton power flow
%   SUCCESS = SOLVE_OPF(OBJ, MPC)
%
%   Inputs
%       OBJ - 
%       MPC - 
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

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(mpc.bus, mpc.gen);

%% create opf_model object
om = opf_model(mpc);
nb = size(mpc.bus, 1);
ng = size(mpc.gen, 1);
baseMVA = mpc.baseMVA;

%% set up variables
vars = horzcat(obj.model_vvars(), obj.model_zvars());
for vtype = vars
    [v0, vl, vu] = obj.params_var(vtype{1});
    om.add_var(vtype{1}, length(v0), v0, vl, vu);
end

%% power balance constraints
fcn_mis = @(x)opf_power_balance_fcn(obj, x);
hess_mis = @(x, lam)opf_power_balance_hess(obj, x, lam);
om.add_nln_constraint({'Pmis', 'Qmis'}, [nb;nb], 1, fcn_mis, hess_mis, vars);

%% find/prepare polynomial generator costs
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
cpg = [];
cqg = [];
[pcost qcost] = pqcost(mpc.gencost, ng);
ip0 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 1);   %% constant
ip1 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 2);   %% linear
ip2 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 3);   %% quadratic
ip3 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) > 3);    %% cubic or greater
if ~isempty(ip2) || ~isempty(ip1) || ~isempty(ip0)
    kpg = zeros(ng, 1);
    cpg = zeros(ng, 1);
    if ~isempty(ip2)
        Qpg = zeros(ng, 1);
        Qpg(ip2) = 2 * pcost(ip2, COST) * baseMVA^2;
        cpg(ip2) = cpg(ip2) + pcost(ip2, COST+1) * baseMVA;
        kpg(ip2) = kpg(ip2) + pcost(ip2, COST+2);
    else
        Qpg = [];   %% no quadratic terms
    end
    if ~isempty(ip1)
        cpg(ip1) = cpg(ip1) + pcost(ip1, COST) * baseMVA;
        kpg(ip1) = kpg(ip1) + pcost(ip1, COST+1);
    end
    if ~isempty(ip0)
        kpg(ip0) = kpg(ip0) + pcost(ip0, COST);
    end
end
nz_extra = om.getN('var', 'zr') - ng;
if nz_extra
    cpg = [cpg; zeros(nz_extra, 1)];
    kpg = [kpg; zeros(nz_extra, 1)];
end
om.add_quad_cost('polPg', Qpg, cpg, kpg, {'zr'});

% from mipsopf_solver()
%
%% options
opt = mpopt.mips;
opt.verbose = mpopt.verbose;
if opt.feastol == 0
    opt.feastol = mpopt.opf.violation;  %% = MPOPT.opf.violation by default
end
if ~isfield(opt, 'cost_mult') || isempty(opt.cost_mult)
    opt.cost_mult = 1e-4;
end

%% unpack data
mpc = om.get_mpc();
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
[vv, ll, nne, nni] = om.get_idx();

%% problem dimensions
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ny = om.getN('var', 'y');   %% number of piece-wise linear costs

%% linear constraints
[A, l, u] = om.params_lin_constraint();

%% bounds on optimization vars
[x0, xmin, xmax] = om.params_var();

%%-----  run opf  -----
f_fcn = @(x)opf_costfcn(x, om);
gh_fcn = @(x)opf_consfcn(x, om);
hess_fcn = @(x, lambda, cost_mult)opf_hessfcn(x, lambda, cost_mult, om);
[x, f, info, Output, Lambda] = ...
  mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt);
success = (info > 0);

i = Output.iterations;
