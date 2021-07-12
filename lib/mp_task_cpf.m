classdef mp_task_cpf < mp_task_pf
%MP_TASK_CPF  MATPOWER task for continuation power flow (CPF).
%   MP_TASK_CPF provides implementation for continuation power flow problem.
%
%   Properties
%       warmstart
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        warmstart   %% warm start data
    end

    methods
        %% constructor
        function obj = mp_task_cpf()
            %% call parent constructor
            obj@mp_task_pf();

            obj.tag = 'CPF';
            obj.name = 'Continuation Power Flow';
        end

        %%-----  task methods  -----
        function [d, mpopt] = run_pre(obj, d, mpopt)
            if ~isa(d, 'mp_data')
                if ~iscell(d) || length(d) < 2
                    error('mp_task_cpf/run_pre: input cases must be provided in a 2-element cell array, specifying the base and target cases, respectively')
                end
                d{1} = run_pre@mp_task_pf(obj, d{1}, mpopt);
                d{2} = run_pre@mp_task_pf(obj, d{2}, mpopt);
            end
        end

        function [results, success] = legacy_post_run(obj, mpopt)
            success = obj.success;
            results = obj.dm.mpc;
            if obj.nm.np ~= 0
                if obj.dm.elements.is_index_name('shunt')
                    results = obj.dmc.elements.shunt.export( ...
                        obj.dm.elements.shunt, results);
                end
                results = obj.dmc.elements.load.export( ...
                    obj.dm.elements.load, results);
                results = obj.dmc.elements.branch.export( ...
                    obj.dm.elements.branch, results, ...
                    {'ql_to', 'pl_to', 'ql_fr', 'pl_fr'});
                results = obj.dmc.elements.gen.export( ...
                    obj.dm.elements.gen, results, ...
                    {'pg', 'qg'});
%                 results = obj.dmc.export(obj.dm, obj.dm.mpc);
            end
        end

        function [mm, nm, dm] = next_mm(obj, mm, nm, dm, mpopt)
            %% return new math model, or empty matrix if finished
            if isfield(mm.soln.output, 'warmstart')
                %% get warmstart info
                ad = mm.get_userdata('aux_data');
                ws = mm.soln.output.warmstart;

                %% save parameter lambda and solved voltages
                %% for current & prev step
                ws.clam = ws.x(end);
                ws.plam = ws.xp(end);
                [ws.cV, ~] = nm.cpf_convert_x(ws.x, ad);
                [ws.pV, ~] = nm.cpf_convert_x(ws.xp, ad);

                %% expand tangent z to all nodes + lambda, for cur & prev step
                [ws.z, ws.zp] = nm.cpf_expand_z_warmstart(ad, ws.z, ws.zp);

                %% set warmstart for next math model
                obj.warmstart = ws;

                %% save updated target models
                nm.userdata.target = ws.nmt;
                dm.userdata.target = ws.dmt;

                %% update network model with current solution
                obj.nm = obj.network_model_update(mm, nm);

                %% update data model voltages only
                %% preserve original base/target specifications
                for k = 1:obj.nm.node.NS
                    nme = obj.nm.elements.(obj.nm.node.order(k).name);
                    nme.pf_data_model_update(mm, obj.nm, obj.dm, mpopt);
                end

                %% create new math model
                mm = obj.math_model_build(nm, dm, mpopt);
            else
                mm = [];
            end
        end

        %%-----  data model converter methods  -----
        function dmc_class = dm_converter_class(obj, d, mpopt)
            if iscell(d) && length(d) == 2
                dmc_class = dm_converter_class@mp_task(obj, d{1}, mpopt);
            else
                error('mp_task_cpf/dm_converter_class: d must be 2-element cell array');
            end
        end

        %%-----  data model methods  -----
        function dm = data_model_build(obj, d, dmc, mpopt)
            if iscell(d) && length(d) == 2
                dm  = data_model_build@mp_task_pf(obj, d{1}, dmc, mpopt);
                dmt = data_model_build@mp_task_pf(obj, d{2}, dmc, mpopt);
                dm.userdata.target = dmt;
            else
                error('mp_task_cpf/data_model_build: d must be 2-element cell array');
            end
        end

        %%-----  network model methods  -----
        function nm = network_model_build(obj, dm, mpopt)
            dmt = dm.userdata.target;
            nm  = network_model_build@mp_task_pf(obj, dm,  mpopt);
            nmt = network_model_build@mp_task_pf(obj, dmt, mpopt);
            nm.userdata.target = nmt;
        end

        %%-----  mathematical model methods  -----
        function mm_class = math_model_class(obj, nm, dm, mpopt)
            mm_class = @mp_math_cpf;
        end

        function opt = math_model_opt(obj, mm, nm, dm, mpopt)
            opt = math_model_opt@mp_task(obj, mm, nm, dm, mpopt);

            %% add the warmstart options, if available
            if ~isempty(obj.warmstart)
                opt = nm.cpf_solve_opts_warmstart(opt, obj.warmstart, mm);
                obj.warmstart = [];     %% delete warmstart data from task
            end
        end
    end     %% methods
end         %% classdef
