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

f = {'PS', 'PI', 'CS', 'CI'};
formulation = [0 0; 0 1; 1 0; 1 1];
cases = {'t_case9_gizmo', 'case4gs', 'case4_dist', 'case9', 'case14', 'case57', 'case300'};
% cases = {'t_case9_gizmo', 'case14'};
% cases = {'case4gs', 't_case9_gizmo', 'case6ww'};
% cases = {'case4_dist'};
% cases = {'case3r'};
% cases = {'case5'};
% cases = {'case9'};
% cases = {'case14'};
% cases = {'t_case9_test'};
% cases = {'t_case9_gizmo'};
% cases = {'case4gs'};

t_begin(6*length(cases)*length(f), quiet);

define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

mpopt = mpoption('out.all', 0, 'verbose', 0, 'pf.tol', 1e-10);
mpopt0 = mpopt;
mpopt0.exp.mpe = 0;
mpopt.exp.data_model_class = @mp_data_mpc2_node_test;
mpopt.exp.network_model_class = @mp_network_acps_node_test;

for k = 1:length(cases)
    for j = 1:4
        mpopt = mpoption(mpopt, 'pf.v_cartesian', formulation(j, 1), ...
                                'pf.current_balance', formulation(j, 2));
        mpopt0 = mpoption(mpopt0, 'pf.v_cartesian', formulation(j, 1), ...
                                'pf.current_balance', formulation(j, 2));
        t = sprintf('PF (%s) - %s - ', f{j}, cases{k});
        mpc = loadcase(cases{k});
        if k == 1
            mpc.bus(2, BS) = 1;
        end

        r = runpf(mpc, mpopt0);
        eVm = r.bus(:, VM);
        eVa = r.bus(:, VA);
        ePg = r.gen(:, PG);
        eQg = r.gen(:, QG);
%         pf0 = run_pf(mpc, mpopt0);
%         eVa = pf0.dm.mpc.bus(:, VA);
%         eVm = pf0.dm.mpc.bus(:, VM);
%         ePg = pf0.dm.mpc.gen(:, PG);
%         eQg = pf0.dm.mpc.gen(:, QG);
        t_ok(r.success, [t 'success 1']);

        pf = run_pf(mpc, mpopt);
        Va = pf.dm.mpc.bus(:, VA);
        Vm = pf.dm.mpc.bus(:, VM);
        Pg = pf.dm.mpc.gen(:, PG);
        Qg = pf.dm.mpc.gen(:, QG);
        t_ok(pf.mm.soln.eflag, [t 'success 2']);
        t_is(Va, eVa, 9, [t 'Va']);
        t_is(Vm, eVm, 9, [t 'Vm']);
        t_is(Pg, ePg, 9, [t 'Pg']);
        t_is(Qg, eQg, 9, [t 'Qg']);

%         t = sprintf('OPF (%s) - %s - ', f{j}, cases{k});
%         r = runopf(mpc, mpopt0);
%         t_ok(r.success, [t 'success 1']);
%         eVm = r.bus(:, VM);
%         eVa = r.bus(:, VA);
%         ePg = r.gen(:, PG);
%         eQg = r.gen(:, QG);
%     
%         opf = run_opf(mpc, mpopt);
%         Va = opf.dm.mpc.bus(:, VA);
%         Vm = opf.dm.mpc.bus(:, VM);
%         Pg = opf.dm.mpc.gen(:, PG);
%         Qg = opf.dm.mpc.gen(:, QG);
%         t_ok(opf.mm.soln.eflag, [t 'success 2']);
%         t_is(Va, eVa, 10, [t 'Va']);
%         t_is(Vm, eVm, 10, [t 'Vm']);
%         t_is(Pg, ePg, 10, [t 'Pg']);
%         t_is(Qg, eQg, 10, [t 'Qg']);
    end
end

t_end;

if nargout
    pf_ = pf;
%     if nargout > 1
%         opf_ = opf;
%     end
end
