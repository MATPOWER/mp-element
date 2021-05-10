function obj = t_node_test(quiet)
%T_NODE_TEST  Tests for network model with multipe node-creating elements.

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

t_begin(3, quiet);

define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

casefile = 't_case9_gizmo';
mpc = loadcase(casefile);
mpc.bus(2, BS) = 1;

mpopt = mpoption('out.all', -1, 'verbose', 2);
mpopt.exp.data_model_class = @mp_data_mpc2_node_test;
mpopt.exp.network_model_class = @mp_network_acps_node_test;

% r = runpf(mpc, mpopt);
% t_ok(r.success, 'success');

pf = run_pf(mpc, mpopt);

Va = pf.dm.mpc.bus(:, VA);
Vm = pf.dm.mpc.bus(:, VM);
Pg = pf.dm.mpc.gen(:, PG);
Qg = pf.dm.mpc.gen(:, QG);
V = [Vm Va];
Sg = [Pg Qg];
eV = [
   1.000000000000000                   0;
   1.000000000000000   9.668741126628102;
   1.000000000000000   4.771073237177304;
   0.987006852391906  -2.406643919519414;
   0.975472177085053  -4.017264326707555;
   1.003375436452801   1.925601686828552;
   0.985644881724947   0.621544555388917;
   0.996185245809070   3.799120192692305;
   0.957621040429904  -4.349933576561018;
];
eSg = 100*[
   0.719547015892220   0.240689577727586;
   0.850000000000000  -0.036490255342106;
   1.630000000000000   0.134601195311253;
];
t_ok(pf.mm.soln.eflag, 'success');
t_is(V, eV, 6, 'V');
t_is(Sg, eSg, 6, 'Sg');

t_end;

if nargout
    obj = pf;
end
