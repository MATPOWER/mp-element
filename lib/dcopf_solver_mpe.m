function [results, success, raw] = dcopf_solver_mpe(opf, mpopt)

%% unpack data
mm = opf.mm;
mpc = opf.dm.mpc;

%% problem dimensions
ny = mm.getN('var', 'y');   %% number of piece-wise linear costs

lambda = mm.soln.lambda;
mu = struct( ...
  'var', struct('l', lambda.lower, 'u', lambda.upper), ...
  'lin', struct('l', lambda.mu_l, 'u', lambda.mu_u) );

results = mpc;
[results.om, results.x, results.mu, results.f] = deal(mm, mm.soln.x, mu, mm.soln.f);

pimul = [ ...
  results.mu.lin.l - results.mu.lin.u;
 -ones(ny>0, 1);    %% dummy entry corresponding to linear cost row in A (in MINOS)
  results.mu.var.l - results.mu.var.u;
];
success = opf.success;
raw = struct('xr', mm.soln.x, 'pimul', pimul, 'info', mm.soln.eflag, 'output', mm.soln.output);
