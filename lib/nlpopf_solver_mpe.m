function [results, success, raw] = nlpopf_solver_mpe(opf, mpopt)

[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% unpack data
mm = opf.mm;
mpc = opf.dm.mpc;

%% problem dimensions
nb = size(mpc.bus, 1);      %% number of buses
nl = size(mpc.branch, 1);   %% number of branches
ny = mm.getN('var', 'y');   %% number of piece-wise linear costs

muSf = mpc.branch(:, MU_SF) * mpc.baseMVA;
muSt = mpc.branch(:, MU_ST) * mpc.baseMVA;

%% package up results
nlnN = 2*nb + 2*nl;     %% because muSf and muSt are nl x 1, not nl2 x 1

%% extract multipliers for nonlinear constraints
lambda = mm.soln.lambda;
kl = find(lambda.eqnonlin(1:2*nb) < 0);
ku = find(lambda.eqnonlin(1:2*nb) > 0);
nl_mu_l = zeros(nlnN, 1);
nl_mu_u = [zeros(2*nb, 1); muSf; muSt];
nl_mu_l(kl) = -lambda.eqnonlin(kl);
nl_mu_u(ku) =  lambda.eqnonlin(ku);

if isfield(lambda, 'ineqnonlin')
    lam_nli = lambda.ineqnonlin;
else
    lam_nli = [];
end

mu = struct( ...
  'var', struct('l', lambda.lower, 'u', lambda.upper), ...
  'nln', struct('l', nl_mu_l, 'u', nl_mu_u), ...
  'nle', lambda.eqnonlin, ...
  'nli', lam_nli, ...
  'lin', struct('l', lambda.mu_l, 'u', lambda.mu_u) );

results = mpc;
[results.om, results.x, results.mu, results.f] = deal(mm, mm.soln.x, mu, mm.soln.f);

pimul = [ ...
  results.mu.nln.l - results.mu.nln.u;
  results.mu.lin.l - results.mu.lin.u;
 -ones(ny>0, 1);    %% dummy entry corresponding to linear cost row in A (in MINOS)
  results.mu.var.l - results.mu.var.u;
];
success = opf.success;
raw = struct('xr', mm.soln.x, 'pimul', pimul, 'info', mm.soln.eflag, 'output', mm.soln.output);
