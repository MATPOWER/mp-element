classdef mp_task_cpf < mp_task_pf
%MP_TASK_CPF  MATPOWER task for continuation power flow (CPF).
%   MP_TASK_CPF provides implementation for continuation power flow problem.
%
%   Properties
%       dmt
%       nmt
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
        dmt         %% data model for target case
        nmt         %% network model for target case
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

        function [mm, nm, dm] = next_mm(obj, mm, nm, dm, mpopt)
            %% return new math model, or empty matrix if finished
            if isfield(mm.soln.output, 'warmstart')
                %% get warmstart info
                ad = mm.get_userdata('aux_data');
                ws = mm.soln.output.warmstart;

                %% save solved voltages and lambda for current & prev step
                ws.clam = ws.x(end);
                ws.plam = ws.xp(end);
                [ws.cV, ~] = nm.cpf_convert_x(ws.x, ad);
                [ws.pV, ~] = nm.cpf_convert_x(ws.xp, ad);

                %% expand tangent z to all nodes + lam, for current & prev step
                z  = ws.z;
                zp = ws.zp;
                ws.z  = zeros(nm.nv, 1);
                ws.zp = zeros(nm.nv, 1);
                k = [ad.pv; ad.pq; nm.nv/2 + ad.pq; nm.nv+1];
                ws.z(k)  = z;
                ws.zp(k) = zp;
                obj.warmstart = ws;

                %% save updated target models
                obj.nmt = ws.nmt;
                obj.dmt = ws.dmt;

                %% update network model with current solution
                obj.nm = obj.network_model_update(mm, nm);

                %% update data model voltages only
                %% preserve original base/target specifications
                for k = 1:obj.nm.node.NS
                    nme = obj.nm.elm_by_name(obj.nm.node.order(k).name);
                    nme.pf_data_model_update(mm, obj.nm, obj.dm, mpopt);
                end

                %% create new math model
                mm = obj.math_model_build(nm, dm, mpopt);
            else
                mm = [];
            end
        end

        %%-----  data model methods  -----
        function dm = data_model_build(obj, d, mpopt)
            if iscell(d) && length(d) == 2
                dm      = data_model_build@mp_task_pf(obj, d{1}, mpopt);
                obj.dmt = data_model_build@mp_task_pf(obj, d{2}, mpopt);
            else
                error('mp_task_cpf: data_model_build: d must be 2-element cell array');
            end
        end

        function dm = data_model_update(obj, mm, nm, dm, mpopt)
            nm.cpf_data_model_update(mm, nm, dm, mpopt);
        end


        %%-----  network model methods  -----
        function nm = network_model_build(obj, dm, mpopt)
            nm      = network_model_build@mp_task_pf(obj, dm, mpopt);
            obj.nmt = network_model_build@mp_task_pf(obj, obj.dmt, mpopt);
        end

        function nm = network_model_convert_x(obj, mm, nm)
            %% convert solved state from math model to network model soln
            [nm.soln.v, nm.soln.z, nm.soln.x] = nm.cpf_convert_x(mm.soln.x, ...
                                                mm.get_userdata('aux_data'));
        end


        %%-----  mathematical model methods  -----
        function mm = math_model_build_pre(obj, mm, nm, dm, mpopt)
            mm.userdata.aux_data = nm.cpf_aux_data(obj.nmt, dm, obj.dmt, mpopt);
        end

        function obj = math_model_add_vars(obj, mm, nm, dm, mpopt)
            nm.cpf_add_vars(mm, nm, dm, mpopt);
        end

        function obj = math_model_add_constraints(obj, mm, nm, dm, mpopt)
            nm.cpf_add_constraints(mm, nm, dm, mpopt);
        end

        function opt = math_model_opt(obj, mm, nm, dm, mpopt)
            opt = nm.cpf_solve_opts(mm, dm, mpopt);

            %% add the warmstart options, if available
            if ~isempty(obj.warmstart)
                opt = nm.cpf_solve_opts_warmstart(opt, obj.warmstart, mm);
                obj.warmstart = [];     %% delete warmstart data from task
            end
            obj.mm_opt = opt;
        end
    end     %% methods
end         %% classdef
