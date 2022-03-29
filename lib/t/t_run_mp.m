function tssk = t_run_mp(quiet)
%T_RUN_MP  Tests for RUN_MP and simple creation and solve of models.

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
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

t_begin(25, quiet);

casefile = 'case9';
casefilet = 'case9target';
opt = struct('verbose', 0);
mpopt = mpoption(opt);
mpopt = mpoption(mpopt, 'out.all', 0);

tasks = {'PF', 'CPF', 'OPF'};
task_classes = {@mp_task_pf, @mp_task_cpf, @mp_task_opf};
d = {casefile, {casefile, casefilet}, casefile};
dmc_class = @mp_dm_converter_mpc2;
dm_classes = {
    @mp_data,
    @mp_data_cpf,
    @mp_data_opf
};
nm_class = @mp_network_acp;
mm_classes = {
    @mp_math_pf_acps,
    @mp_math_cpf_acps,
    @mp_math_opf_acps
};

t = 'build data model converter';
dmc = dmc_class().build();
t_ok(isa(dmc, 'mp_dm_converter'), t);

for k = 1:length(tasks)
    t = sprintf('%s : ', tasks{k});

    %% build data model
    dm_class = dm_classes{k};
    dm = dm_class().build('case9', dmc);
    t_ok(isa(dm, 'mp_data'), [t 'build data model']);

    %% build network model
    nm = nm_class().build(dm);
    t_ok(isa(nm, 'mp_network'), [t 'build network model']);

    if strcmp(tasks{k}, 'CPF')
        dm.userdata.target = dm_class().build('case9target', dmc);
        nm.userdata.target = nm_class().build(dm.userdata.target);
    end

    %% build math model
    mm_class = mm_classes{k};
    mm = mm_class().build(nm, dm, mpopt);
    t_ok(isa(mm, 'mp_math'), [t 'build math model']);

    %% solve math model
    % opt = mm.solve_opts(nm, dm, mpopt);
    mm.solve(opt);
    t_is(mm.soln.eflag, 1, 12, [t 'solve math model']);

    %% network model solution
    nm = mm.network_model_x_soln(nm);
    t_ok(isfield(nm.soln, 'x'), [t 'network model x soln']);
    nm.port_inj_soln();
    t_ok(isfield(nm.soln, 'gs_'), [t 'network model port inj soln']);

    %% data model update
    dm = mm.data_model_update(nm, dm, mpopt);
    t_ok(all(dm.elements.branch.tab.pl_fr ~= 0), [t 'data model update']);

    %% run_mp
    tsk = run_mp(task_classes{k}, d{k}, mpopt);
    t_is(tsk.success, 1, 12, [t 'run_mp : success']);
end

if nargout  %% set output arg
    tssk = tsk;
end

t_end;
