function tssk = t_run_mp(quiet)
%T_RUN_MP  Tests for RUN_MP and simple creation and solve of models.

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

% define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

t_begin(21, quiet);

casefile = 'case9';
casefilet = 'case9target';
opt = struct('verbose', 0);
mpopt = mpoption(opt);
mpopt = mpoption(mpopt, 'out.all', 0);

dmc_class = @mp_dm_converter_mpc2;
dm_class = @mp_data_mpc2;
nm_class = @mp_network_acps;
mm_classes = {
    @mp_math_pf,
    @mp_math_cpf,
    @mp_math_opf
};
tasks = {'PF', 'CPF', 'OPF'};
d = {casefile, {casefile, casefilet}, casefile};

t = 'build data model converter';
dmc = dmc_class().build();
t_ok(isa(dmc, 'mp_dm_converter'), t);

t = 'build data model';
dm = dm_class().build('case9', dmc);
t_ok(isa(dm, 'mp_data'), t);

t = 'build network model';
nm = nm_class().build(dm);
t_ok(isa(nm, 'mp_network'), t);

dm.userdata.target = dm_class().build('case9target', dmc);
nm.userdata.target = nm_class().build(dm.userdata.target);

for k = 1:length(mm_classes)
    t = sprintf('%s : ', tasks{k});
    mm_class = mm_classes{k};
    mm = mm_class().build(nm, dm, mpopt);
    t_ok(isa(mm, 'mp_math'), [t 'build math model']);

    % opt = mm.solve_opts(nm, dm, mpopt);
    mm.solve(opt);
    t_is(mm.soln.eflag, 1, 12, [t 'solve']);
    nm = mm.network_model_x_soln(nm);
    t_ok(isfield(nm.soln, 'x'), [t 'network model soln']);
    nm.port_inj_soln();
% nm.soln
    t_ok(isfield(nm.soln, 'gs_'), [t 'network model port inj soln']);
% keyboard
    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
        TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
        ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
% Pf = dm.mpc.branch(:, PF)
    dm = mm.data_model_update(nm, dm, mpopt);
% Pf = dm.mpc.branch(:, PF)
    t_ok(all(dm.mpc.branch(:, PF) ~= 0), [t 'data model update']);

    tsk = run_mp(tasks{k}, d{k}, mpopt);
    t_is(tsk.success, 1, 12, [t 'run_mp : success']);
end

if nargout  %% set output arg
    tssk = tsk;
end

t_end;
