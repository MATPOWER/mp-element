function [pf_, opf_] = t_node_test(quiet)
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

cases = {'t_case9_gizmo', 'case4gs', 'case4_dist', 'case9', 'case14', 'case57', 'case300'};

t_begin(12*length(cases), quiet);

define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

mpopt = mpoption('out.all', 0, 'verbose', 0, 'pf.tol', 1e-10);
mpopt = mpoption(mpopt, 'opf.ignore_angle_lim', 1);
mpopt0 = mpopt;
mpopt0.exp.mpe = 0;
mpopt.exp.dm_converter_class = @mp_dm_converter_mpc2_node_test;
mpopt.exp.data_model_class = @mp_data_mpc2_node_test;
mpopt.exp.network_model_class = @mp_network_acps_node_test;

for k = 1:length(cases)
    t = sprintf('PF - %s - ', cases{k});
    mpc = loadcase(cases{k});
    if strcmp(cases{k}, 't_case9_gizmo')
        mpc.bus(2, BS) = 1;
    end
    if ~isfield(mpc, 'gencost')
        ng = size(mpc.gen, 1);
        mpc.gencost = ones(ng, 1) * [2 0 0 2 1 0];
        if strcmp(cases{k}, 'case4gs')
            mpc.gen(:, PMAX) = 320;
            mpc.gen(:, QMAX) = 200;
        end
    end

    have_feature('mp_element', 0);
    r = runpf(mpc, mpopt0);
    have_feature('mp_element', 1);
    eVm = r.bus(:, VM);
    eVa = r.bus(:, VA);
    ePg = r.gen(:, PG);
    eQg = r.gen(:, QG);
    t_ok(r.success, [t 'success 1']);

    pf = run_pf(mpc, mpopt);
    if pf.nm.np ~= 0
        pf.dm.userdata.mpc = pf.dmc.export(pf.dm, pf.dm.userdata.mpc, pf.tag);
    end
    Va = pf.dm.userdata.mpc.bus(:, VA);
    Vm = pf.dm.userdata.mpc.bus(:, VM);
    Pg = pf.dm.userdata.mpc.gen(:, PG);
    Qg = pf.dm.userdata.mpc.gen(:, QG);
    t_ok(pf.success, [t 'success 2']);
    t_is(Va, eVa, 9, [t 'Va']);
    t_is(Vm, eVm, 9, [t 'Vm']);
    t_is(Pg, ePg, 9, [t 'Pg']);
    t_is(Qg, eQg, 9, [t 'Qg']);

    t = sprintf('OPF - %s - ', cases{k});
    have_feature('mp_element', 0);
    r = runopf(mpc, mpopt0);
    have_feature('mp_element', 1);
    t_ok(r.success, [t 'success 1']);
    eVm = r.bus(:, VM);
    eVa = r.bus(:, VA);
    ePg = r.gen(:, PG);
    eQg = r.gen(:, QG);

    opf = run_opf(mpc, mpopt);
    if opf.nm.np ~= 0
        opf.dm.userdata.mpc = opf.dmc.export(opf.dm, opf.dm.userdata.mpc, opf.tag);
    end
    Va = opf.dm.userdata.mpc.bus(:, VA);
    Vm = opf.dm.userdata.mpc.bus(:, VM);
    Pg = opf.dm.userdata.mpc.gen(:, PG);
    Qg = opf.dm.userdata.mpc.gen(:, QG);
    t_ok(opf.success, [t 'success 2']);
    t_is(Va, eVa, 9, [t 'Va']);
    t_is(Vm, eVm, 9, [t 'Vm']);
    t_is(Pg, ePg, 9, [t 'Pg']);
    t_is(Qg, eQg, 9, [t 'Qg']);
end

t_end;

if nargout
    pf_ = pf;
    if nargout > 1
        opf_ = opf;
    end
end
