classdef task_opf_legacy < mp.task_opf & mp.task_shared_legacy
%MP.TASK_OPF_LEGACY  MATPOWER task for legacy optimal power flow (OPF).
%   MP.TASK_OPF_LEGACY provides implementation for optimal power flow problem.
%
%   Properties
%       ?
%
%   Methods
%       ?
%
%   See also MP.TASK_OPF

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
    end

    methods
        %%-----  task methods  -----
        function [d, mpopt] = run_pre(obj, d, mpopt)
            [d, mpopt] = obj.run_pre_legacy(d, mpopt);
            [d, mpopt] = run_pre@mp.task_opf(obj, d, mpopt);
        end

        function obj = run_post(obj, mm, nm, dm, mpopt);
            if obj.nm.np ~= 0
                obj.dm.source = obj.dmc.export(obj.dm, obj.dm.source);
            end
        end

        %%-----  data model converter methods  -----
        function dmc_class = dm_converter_class_mpc2_default(obj)
            dmc_class = @mp.dm_converter_mpc2_legacy;
        end

        %%-----  data model methods  -----
        function dm = data_model_build_post(obj, dm, dmc, mpopt)
            %% call parent
            dm = data_model_build_post@mp.task_opf(obj, dm, dmc, mpopt);

            %% pre-process inputs for legacy user vars, constraints, costs
            dm = dmc.legacy_user_mod_inputs(dm, mpopt, obj.dc);
        end

        %%-----  mathematical model methods  -----
        function mm_class = math_model_class_default(obj, nm, dm, mpopt)
            %% mp.mm_shared_opf_legacy (compatible with opf_model) is required
            %% to support legacy cost functions and callback functions that
            %% expect to find mpc in mm.mpc.
            switch upper(mpopt.model)
                case 'AC'
                    if mpopt.opf.v_cartesian
                        if mpopt.opf.current_balance
                            mm_class = @mp.math_model_opf_acci_legacy;
                        else
                            mm_class = @mp.math_model_opf_accs_legacy;
                        end
                    else
                        if mpopt.opf.current_balance
                            mm_class = @mp.math_model_opf_acpi_legacy;
                        else
                            mm_class = @mp.math_model_opf_acps_legacy;
                        end
                    end
                case 'DC'
                    mm_class = @mp.math_model_opf_dc_legacy;
            end
        end

        function [results, success, raw] = legacy_post_run(obj, mpopt)
            %% from opf_execute(), dcopf_solver(), nlpopf_solver()

            %% unpack data
            mm = obj.mm;
            success = obj.success;
            results = obj.dm.source;
            results.success = success;

            if mm.getN('var')
                lambda = mm.soln.lambda;
                dc  = obj.dc;
                ny = mm.getN('var', 'y');   %% number of piece-wise linear costs

                %% get indexing
                [vv, ll, nne, nni] = mm.get_idx();

                if dc
                    mu = struct( ...
                      'var', struct('l', lambda.lower, 'u', lambda.upper), ...
                      'lin', struct('l', lambda.mu_l, 'u', lambda.mu_u) );

                    pimul = [ ...
                      mu.lin.l - mu.lin.u;
                     -ones(ny>0, 1);    %% dummy entry corresponding to linear cost row in A (in MINOS)
                      mu.var.l - mu.var.u;
                    ];
                else
                    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                        TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                        ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

                    %% problem dimensions
                    nb = size(results.bus, 1);      %% number of buses
                    nl = size(results.branch, 1);   %% number of branches

                    muSf = results.branch(:, MU_SF) * results.baseMVA;
                    muSt = results.branch(:, MU_ST) * results.baseMVA;

                    %% package up results
                    nlnN = 2*nb + 2*nl;     %% because muSf and muSt are nl x 1, not nl2 x 1

                    %% extract multipliers for nonlinear constraints
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

                    pimul = [ ...
                        mu.nln.l - mu.nln.u;
                        mu.lin.l - mu.lin.u;
                       -ones(ny>0, 1);    %% dummy entry corresponding to linear cost row in A (in MINOS)
                        mu.var.l - mu.var.u;
                    ];
                end %% if dc

                [results.om, results.x, results.mu, results.f] = ...
                    deal(mm, mm.soln.x, mu, mm.soln.f);
                raw = struct('xr', mm.soln.x, 'pimul', pimul, 'info', mm.soln.eflag, ...
                            'output', mm.soln.output, 'task', obj);

                if ~isfield(raw, 'output') || ~isfield(raw.output, 'alg') || isempty(raw.output.alg)
                    raw.output.alg = upper(mpopt.opf.ac.solver);
                end

                if success
                    if ~dc
                        %% compute g, dg, f, df, d2f if requested by opf.return_raw_der = 1
                        if mpopt.opf.return_raw_der
                            %% move from results to raw if using v4.0 of MINOPF or TSPOPF
                            if isfield(results, 'dg')
                                raw.dg = results.dg;
                                raw.g = results.g;
                            end
                            %% compute g, dg, unless already done by post-v4.0 MINOPF or TSPOPF
                            if ~isfield(raw, 'dg')
                                [g, geq, dg, dgeq] = nlp_consfcn(mm, results.x);
                                raw.g = [ geq; g];
                                raw.dg = [ dgeq'; dg'];   %% true Jacobian organization
                            end
                            %% compute df, d2f
                            [f, df, d2f] = nlp_costfcn(mm, results.x);
                            raw.df = df;
                            raw.d2f = d2f;
                        end
                    end

                    %% delete g and dg fields from results if using v4.0 of MINOPF or TSPOPF
                    if isfield(results, 'dg')
                        rmfield(results, 'dg');
                        rmfield(results, 'g');
                    end
                else
                    raw.output.message = obj.message;

                    %% assign empty g, dg, f, df, d2f if requested by opf.return_raw_der = 1
                    if ~dc && mpopt.opf.return_raw_der
                        raw.dg = [];
                        raw.g = [];
                        raw.df = [];
                        raw.d2f = [];
                    end
                end

                %% assign values and limit shadow prices for variables
                om_var_order = mm.get('var', 'order');
                for k = 1:length(om_var_order)
                    name = om_var_order(k).name;
                    if mm.getN('var', name)
                        idx = vv.i1.(name):vv.iN.(name);
                        results.var.val.(name) = results.x(idx);
                        results.var.mu.l.(name) = results.mu.var.l(idx);
                        results.var.mu.u.(name) = results.mu.var.u(idx);
                    end
                end

                %% assign shadow prices for linear constraints
                om_lin_order = mm.get('lin', 'order');
                for k = 1:length(om_lin_order)
                    name = om_lin_order(k).name;
                    if mm.getN('lin', name)
                        idx = ll.i1.(name):ll.iN.(name);
                        results.lin.mu.l.(name) = results.mu.lin.l(idx);
                        results.lin.mu.u.(name) = results.mu.lin.u(idx);
                    end
                end

                %% assign shadow prices for nonlinear constraints
                if ~dc
                    om_nle_order = mm.get('nle', 'order');
                    for k = 1:length(om_nle_order)
                        name = om_nle_order(k).name;
                        if mm.getN('nle', name)
                            results.nle.lambda.(name) = results.mu.nle(nne.i1.(name):nne.iN.(name));
                        end
                    end

                    om_nli_order = mm.get('nli', 'order');
                    for k = 1:length(om_nli_order)
                        name = om_nli_order(k).name;
                        if mm.getN('nli', name)
                            results.nli.mu.(name) = results.mu.nli(nni.i1.(name):nni.iN.(name));
                        end
                    end
                end

                %% assign values for components of quadratic cost
                om_qdc_order = mm.get('qdc', 'order');
                for k = 1:length(om_qdc_order)
                    name = om_qdc_order(k).name;
                    if mm.getN('qdc', name)
                        results.qdc.(name) = mm.eval_quad_cost(results.x, name);
                    end
                end

                %% assign values for components of general nonlinear cost
                om_nlc_order = mm.get('nlc', 'order');
                for k = 1:length(om_nlc_order)
                    name = om_nlc_order(k).name;
                    if mm.getN('nlc', name)
                        results.nlc.(name) = mm.eval_nln_cost(results.x, name);
                    end
                end

                %% assign values for components of legacy user cost
                om_cost_order = mm.get('cost', 'order');
                for k = 1:length(om_cost_order)
                    name = om_cost_order(k).name;
                    if mm.getN('cost', name)
                        results.cost.(name) = mm.eval_legacy_cost(results.x, name);
                    end
                end

                %% if single-block PWL costs were converted to POLY, insert dummy y into x
                %% Note: The "y" portion of x will be nonsense, but everything should at
                %%       least be in the expected locations.
                gen_dmce = obj.dmc.elements.gen;
                if ~isempty(gen_dmce.pwl1)
                    pwl1 = gen_dmce.pwl1;
                else
                    pwl1 = [];
                end
                if ~isempty(pwl1)
                    %% get indexing
                    vv = mm.get_idx();
                    if dc
                        nx = vv.iN.Pg;
                    else
                        nx = vv.iN.Qg;
                    end
                    y = zeros(length(pwl1), 1);
                    raw.xr = [ raw.xr(1:nx); y; raw.xr(nx+1:end)];
                    results.x = [ results.x(1:nx); y; results.x(nx+1:end)];
                end
            else
                raw.output.message = obj.message;
            end
        end
    end     %% methods
end         %% classdef
