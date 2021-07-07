function tssk = t_mp_dm_converter_mpc2(quiet)
%T_MP_DM_CONVERTER_MPC2  Tests for MP_DM_CONVERTER_MPC2.

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

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

casefile = 't_case_ext';
mpc0 = loadcase(casefile);
mpc0.branch(end, 17) = 0;    %% add 4 columns for PF, QF, PT, QT
dmc = mp_dm_converter_mpc2().build();

t_begin(21, quiet);

t = 'dmc constructor : ';
dmc = mp_dm_converter_mpc2();
t_ok(isa(dmc, 'mp_dm_converter'), [t 'class']);
t_is(length(dmc.element_classes), 5, 12, [t '# of element_classes']);
t_ok(isempty(dmc.elements), [t 'elements empty']);

t = 'dmc.build() : ';
dmc.build();
t_is(length(dmc.elements), 5, 12, [t '# of elements']);
t_ok(strcmp(dmc.elements{1}.name, 'bus'), [t 'element{1} is bus']);
t_ok(strcmp(dmc.elements{2}.name, 'gen'), [t 'element{2} is gen']);
t_ok(strcmp(dmc.elements{3}.name, 'load'), [t 'element{3} is load']);
t_ok(strcmp(dmc.elements{4}.name, 'branch'), [t 'element{4} is branch']);
t_ok(strcmp(dmc.elements{5}.name, 'shunt'), [t 'element{5} is shunt']);

t = 'dm constructor : ';
dm = mp_data_mpc2();
t_ok(isa(dm, 'mp_data'), [t 'class']);
t_is(length(dm.element_classes), 5, 12, [t '# of element_classes']);
t_ok(isempty(dm.elements), [t 'elements empty']);

t = 'dm.build(mpc, dmc) : ';
dm.build(mpc0, dmc);
t_is(length(dm.elements), 5, 12, [t '# of elements']);
t_ok(strcmp(dm.elements{1}.name, 'bus'), [t 'element{1} is bus']);
t_ok(strcmp(dm.elements{2}.name, 'gen'), [t 'element{2} is gen']);
t_ok(strcmp(dm.elements{3}.name, 'load'), [t 'element{3} is load']);
t_ok(strcmp(dm.elements{4}.name, 'branch'), [t 'element{4} is branch']);
t_ok(strcmp(dm.elements{5}.name, 'shunt'), [t 'element{5} is shunt']);

t = 'mpc = dmc.export(dm, mpc0) : ';
mpc = dmc.export(dm, mpc0);
t_ok(isequal(mpc, mpc0), [t 'mpc == mpc0']);

t = 'mpc = dmc.export(updated_dm, mpc0) : ';
mpc1 = mpc0;
mpc1.bus([2;3;6;9], VM) = [1.02;1.03;1.06;1.09];
dm.elements.bus.tab.vm([2;3;6;9]) = [1.02;1.03;1.06;1.09];
mpc1.bus([3;5;7], PD) = [30;50;70];
mpc1.bus([3;5;7], QD) = [18;30;42];
load0 = dm.elements.load.tab;
load = [load0; load0(1:2, :)];
k = [4;1;5];
load.uid(k) = [4;1;5];
load.bus(k) = [30;5;6];
load.source_uid(k) = [3;5;7];
load.pd(k) = [30;50;70];
load.qd(k) = [18;30;42];
dm.elements.load.tab = load;
dm.elements.load.rebuild(dm);
mpc = dmc.export(dm, mpc0);
t_ok(~isequal(mpc, mpc0), [t 'mpc ~= mpc0']);
t_ok( isequal(mpc, mpc1), [t 'mpc == updated_mpc']);


% Vm = [mpc0.bus(:, 8) mpc1.bus(:, 8) mpc.bus(:, 8) dm.elements.bus.tab.vm]
% Pd = [mpc0.bus(:, 3) mpc1.bus(:, 3) mpc.bus(:, 3)]
% Pdm = dm.elements.load.tab.pd

% dmc
% dm

% keyboard

% -- -----------   ---------  --------  --------  --------------------
%  1  bus           bus            10         9    dme_bus_mpc2
%  2  gen           gen             4         3    dme_gen_mpc2
%  3  load          bus             3         3    dme_load_mpc2
%  4  branch        branch         10         9    dme_branch_mpc2
%  5  shunt         bus             2         2    dme_shunt


t_end;
