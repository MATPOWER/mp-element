function [results, success, raw] = dcopf_solver_mpe(opf, mpopt)

% %%----- initialization -----
% %% define named indices into data matrices
% [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
%     VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
% [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
%     MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
%     QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
% [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
%     TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
%     ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
% [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% unpack data
mm = opf.mm;
mpc = opf.dm.mpc;
% [baseMVA, bus, gen, branch, gencost] = ...
%     deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
% [vv, ll] = mm.get_idx();

%% problem dimensions
% nb = size(bus, 1);          %% number of buses
% nl = size(branch, 1);       %% number of branches
ny = mm.getN('var', 'y');   %% number of piece-wise linear costs

% [x, f, eflag, output, lambda, success] = deal(mm.soln.x, mm.soln.f, ...
%     mm.soln.eflag, mm.soln.output, mm.soln.lambda, opf.success);
% 
%%-----  calculate return values  -----
% Bf = mm.get_userdata('Bf');
% Pfinj = mm.get_userdata('Pfinj');
% if ~any(isnan(x))
%     %% update solution data
%     Va = x(vv.i1.Va:vv.iN.Va);
%     Pg = x(vv.i1.Pg:vv.iN.Pg);
% 
%     %% update voltages & generator outputs
%     bus(:, VM) = ones(nb, 1);
%     bus(:, VA) = Va * 180/pi;
%     gen(:, PG) = Pg * baseMVA;
% 
%     %% compute branch flows
%     branch(:, [QF, QT]) = zeros(nl, 2);
%     branch(:, PF) = (Bf * Va + Pfinj) * baseMVA;
%     branch(:, PT) = -branch(:, PF);
% end
% 
% %% package up results
% mu_l = lambda.mu_l;
% mu_u = lambda.mu_u;
% muLB = lambda.lower;
% muUB = lambda.upper;
% 
% %% update Lagrange multipliers
% il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
% bus(:, [LAM_P, LAM_Q, MU_VMIN, MU_VMAX]) = zeros(nb, 4);
% gen(:, [MU_PMIN, MU_PMAX, MU_QMIN, MU_QMAX]) = zeros(size(gen, 1), 4);
% branch(:, [MU_SF, MU_ST]) = zeros(nl, 2);
% bus(:, LAM_P)       = (mu_u(ll.i1.Pmis:ll.iN.Pmis) - mu_l(ll.i1.Pmis:ll.iN.Pmis)) / baseMVA;
% if ~isempty(il)
%     branch(il, MU_SF)   = mu_u(ll.i1.Pf:ll.iN.Pf) / baseMVA;
%     branch(il, MU_ST)   = mu_l(ll.i1.Pf:ll.iN.Pf) / baseMVA;
% end
% gen(:, MU_PMIN)     = muLB(vv.i1.Pg:vv.iN.Pg) / baseMVA;
% gen(:, MU_PMAX)     = muUB(vv.i1.Pg:vv.iN.Pg) / baseMVA;

lambda = mm.soln.lambda;
mu = struct( ...
  'var', struct('l', lambda.lower, 'u', lambda.upper), ...
  'lin', struct('l', lambda.mu_l, 'u', lambda.mu_u) );

results = mpc;
% [results.bus, results.branch, results.gen, ...
%     results.om, results.x, results.mu, results.f] = ...
%         deal(bus, branch, gen, mm, x, mu, f);
[results.om, results.x, results.mu, results.f] = deal(mm, mm.soln.x, mu, mm.soln.f);

pimul = [
  results.mu.lin.l - results.mu.lin.u;
 -ones(ny>0, 1);    %% dummy entry corresponding to linear cost row in A (in MINOS)
  results.mu.var.l - results.mu.var.u;
];
success = opf.success;
raw = struct('xr', mm.soln.x, 'pimul', pimul, 'info', mm.soln.eflag, 'output', mm.soln.output);
