classdef mp_task < handle
%MP_TASK  MATPOWER task abstract base class.
%   Each task type (e.g. power flow, CPF, OPF) will inherit from MP_TASK.
%
%   MP_TASK provides properties and methods related to the specific
%   problem specification being solved (e.g. power flow, continuation
%   power flow, optimal power flow, etc.). In particular, it coordinates
%   all interactions between the 3 model layers: data model, network model,
%   and mathematical model.
%
%   Properties
%       dm - data model
%       nm - network model
%       mm - mathematical model
%       mm_opt - solve options for mathematical model
%       tag - task tag - e.g. 'PF', 'CPF', 'OPF'
%       name - task name - e.g. 'Power Flow', etc.
%       i_dm - iteration counter for data model loop
%       i_nm - iteration counter for network model loop
%       i_mm - iteration counter for math model loop
%       success - success flag, 1 - math model solved, 0 - didn't solve
%       message - output message
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        dmc     %% data model converter object
        dm      %% data model object
        nm      %% network model object
        mm      %% mathematical model object
        mm_opt  %% solve options for mathematical model
        tag     %% task tag - e.g. 'PF', 'CPF', 'OPF'
        name    %% task name - e.g. 'Power Flow', etc.
        i_dm    %% iteration counter for data model loop
        i_nm    %% iteration counter for network model loop
        i_mm    %% iteration counter for math model loop
        success %% success flag, 1 - math model solved, 0 - didn't solve
        message %% output message
    end

    methods
        %%-----  task methods  -----
        function obj = run(obj, d, mpopt)
            [d, mpopt] = obj.run_pre(d, mpopt);

            dmc = obj.dm_converter_create(d, mpopt);
            obj.dmc = dmc;

            %% initialize
            obj.i_dm = 0;   %% iteration counter for data model loop
            obj.i_nm = 0;   %% iteration counter for network model loop
            obj.i_mm = 0;   %% iteration counter for math model loop

            %% build data model
            obj.i_dm = obj.i_dm + 1;
            dm = obj.data_model_build(d, dmc, mpopt);
            obj.dm = dm;    %% stash current data model in task object

            while ~isempty(dm)  %% begin data model loop
                %% build network model
                obj.i_nm = obj.i_nm + 1;
                nm = obj.network_model_build(dm, mpopt);
                obj.nm = nm;    %% stash current network model in task object

                while ~isempty(nm)  %% begin network model loop
                    %% build math model
                    obj.i_mm = obj.i_mm + 1;
                    mm = obj.math_model_build(nm, dm, mpopt);
                    obj.mm = mm;    %% stash current math model in task object

                    %% print initial output
                    if mpopt.verbose && obj.i_dm == 1 && obj.i_nm == 1
                        v = mpver('all');
                        fprintf('\nMATPOWER Version %s, %s\n', v.Version, v.Date);
                        fprintf('%s -- %s formulation\n', ...
                            mm.task_name(), mm.form_name());
                    end

                    while ~isempty(mm)  %% begin math model loop
                        if mm.getN('var') == 0  %% model IS empty
                            obj.success = 0;
                            obj.message = sprintf('%s not valid : MATPOWER model contains no connected buses', obj.tag);
                            repeat_mm = 0;
                        else                    %% model is NOT empty
                            %% get solve options
                            mm_opt = obj.math_model_opt(mm, nm, dm, mpopt);
                            obj.mm_opt = mm_opt;    %% stash math model solve
                                                    %% options in task object

                            %% solve mathematical model
                            mm.solve(mm_opt);
                            obj.success = (mm.soln.eflag > 0);
                            if obj.success
                                obj.message = sprintf('%s successful', obj.tag);
                            else
                                obj.message = sprintf('%s failed', obj.tag);
                            end
                        end

                        [mm, nm, dm] = obj.next_mm(mm, nm, dm, mpopt);
                        if ~isempty(mm)
                            obj.mm = mm;
                            obj.i_mm = obj.i_mm + 1;
                        end
                    end                 %% end math model loop
                    mm = obj.mm;        %% use stashed math model below

                    %% update network model with math model solution
                    if nm.np == 0
                        nm = [];
                    else
                        nm = obj.network_model_update(mm, nm);

                        [nm, dm] = obj.next_nm(mm, nm, dm, mpopt);
                        if ~isempty(nm)
                            obj.nm = nm;
                            obj.i_nm = obj.i_nm + 1;
                        end
                    end
                end                 %% end network model loop
                nm = obj.nm;        %% use stashed network model below

                %% update data model with network/math model solution
                dm = mm.data_model_update(nm, dm, mpopt);
                if mpopt.verbose
                    fprintf('%s\n', obj.message);
                end

                dm = obj.next_dm(mm, nm, dm, mpopt);
                if ~isempty(dm)
                    obj.dm = dm;
                    obj.i_dm = obj.i_dm + 1;
                end
            end                 %% end data model loop
            obj.run_post(mm, nm, obj.dm, mpopt);
        end

        function [mm, nm, dm] = next_mm(obj, mm, nm, dm, mpopt)
            %% return new math model, or empty matrix if finished
            mm = [];
        end

        function [nm, dm] = next_nm(obj, mm, nm, dm, mpopt)
            %% return new network model, or empty matrix if finished
            nm = [];
        end

        function dm = next_dm(obj, mm, nm, dm, mpopt)
            %% return new data model, or empty matrix if finished
            dm = [];
        end

        function [d, mpopt] = run_pre(obj, d, mpopt)
        end

        function obj = run_post(obj, mm, nm, dm, mpopt);
        end

        function print_soln(obj, fname, mpopt)
            if mpopt.out.all
                if isempty(fname)
                    fprintf('-- %s print_soln()\n', obj.tag);
                else
                    fprintf('-- %s print_soln(''%s'')\n', obj.tag, fname);
                end
            end
        end

        function save_soln(obj, fname)
            fprintf('-- %s save_soln(''%s'')\n', obj.tag, fname);
        end

        %%-----  data model converter methods  -----
        function dmc_class = dm_converter_class(obj, d, mpopt)
            if isfield(mpopt.exp, 'dm_converter_class') && ...
                    ~isempty(mpopt.exp.dm_converter_class)
                dmc_class = mpopt.exp.dm_converter_class;
            else
                if ismpc2(d)
                    dmc_class = obj.dm_converter_class_mpc2_default();
                else
                    error('mp_task: input data format not recognized');
                end
            end
        end

        function dmc_class = dm_converter_class_mpc2_default(obj)
            dmc_class = @mp_dm_converter_mpc2;
        end

        function dmc = dm_converter_create(obj, d, mpopt)
            if isa(d, 'mp_data')
                dmc = [];
            else
                dmc_class = obj.dm_converter_class(d, mpopt);
                dmc = dmc_class();

                %% add user-supplied elements to dm.element_classes
                if isfield(mpopt.exp, 'dmc_element_classes') && ...
                        ~isempty(mpopt.exp.dmc_element_classes)
                    dmc.modify_element_classes(mpopt.exp.dmc_element_classes);
                end

                dmc.build();
            end
        end

        %%-----  data model methods  -----
        function dm_class = data_model_class(obj, d, mpopt)
            if isfield(mpopt.exp, 'data_model_class') && ...
                    ~isempty(mpopt.exp.data_model_class)
                dm_class = mpopt.exp.data_model_class;
            else
                dm_class = obj.data_model_class_default();
            end
        end

        function dm_class = data_model_class_default(obj)
            dm_class = @mp_data;
        end

        function dm = data_model_create(obj, d, mpopt)
            dm_class = obj.data_model_class(d, mpopt);
            dm = dm_class();
        end

        function dm = data_model_build(obj, d, dmc, mpopt)
            if isa(d, 'mp_data')
                dm = d;
            else
                dm = obj.data_model_create(d, mpopt);
                [dm, d] = obj.data_model_build_pre(dm, d, dmc, mpopt);
                dm.build(d, dmc);
                dm = obj.data_model_build_post(dm, dmc, mpopt);
            end
        end

        function [dm, d] = data_model_build_pre(obj, dm, d, dmc, mpopt)
            %% add user-supplied elements to dm.element_classes
            if isfield(mpopt.exp, 'dm_element_classes') && ...
                    ~isempty(mpopt.exp.dm_element_classes)
                dm.modify_element_classes(mpopt.exp.dm_element_classes);
            end
        end

        function dm = data_model_build_post(obj, dm, dmc, mpopt)
        end

        %%-----  network model methods  -----
        function nm_class = network_model_class(obj, dm, mpopt)
            if isfield(mpopt.exp, 'network_model_class') && ...
                    ~isempty(mpopt.exp.network_model_class)
                nm_class = mpopt.exp.network_model_class;
            else
                nm_class = obj.network_model_class_default(dm, mpopt);
            end
        end

        function nm_class = network_model_class_default(obj, dm, mpopt)
            error('mp_task/network_model_class_default: must be implemented in sub-class');
        end

        function nm = network_model_create(obj, dm, mpopt)
            nm_class = obj.network_model_class(dm, mpopt);
            nm = nm_class().init_set_types();
        end

        function nm = network_model_build(obj, dm, mpopt)
            nm = obj.network_model_create(dm, mpopt);
            nm = obj.network_model_build_pre(nm, dm, mpopt);
            nm.build(dm);
            nm = obj.network_model_build_post(nm, dm, mpopt);
        end
        
        function nm = network_model_build_pre(obj, nm, dm, mpopt)
            %% add user-supplied elements to nm.element_classes
            if isfield(mpopt.exp, 'nm_element_classes') && ...
                    ~isempty(mpopt.exp.nm_element_classes)
                nm.modify_element_classes(mpopt.exp.nm_element_classes);
            end
        end

        function nm = network_model_build_post(obj, nm, dm, mpopt)
        end

        function nm = network_model_x_soln(obj, mm, nm)
            nm = mm.network_model_x_soln(nm);
        end

        function nm = network_model_update(obj, mm, nm)
            %% save network state solution (convert from math model state)
            obj.network_model_x_soln(mm, nm);

            %% save port injection solution
            nm.port_inj_soln();
        end

        %%-----  mathematical model methods  -----
        function mm_class = math_model_class(obj, nm, dm, mpopt)
            if isfield(mpopt.exp, 'math_model_class') && ...
                    ~isempty(mpopt.exp.math_model_class)
                mm_class = mpopt.exp.math_model_class;
            else
                mm_class = obj.math_model_class_default(nm, dm, mpopt);
            end
        end

        function mm = math_model_create(obj, nm, dm, mpopt)
            mm_class = obj.math_model_class(nm, dm, mpopt);
            mm = mm_class().init_set_types();
        end

        function mm = math_model_build(obj, nm, dm, mpopt)
            mm = obj.math_model_create(nm, dm, mpopt);

            if nm.np ~= 0       %% skip for empty model
%                 mm = obj.math_model_build_pre(mm, nm, dm, mpopt);
                mm.build(nm, dm, mpopt);
%                 mm = obj.math_model_build_post(mm, nm, dm, mpopt);
            end
        end

%         function mm = math_model_build_pre(obj, mm, nm, dm, mpopt)
%         end
% 
%         function mm = math_model_build_post(obj, mm, nm, dm, mpopt)
%         end

        function opt = math_model_opt(obj, mm, nm, dm, mpopt)
            opt = mm.solve_opts(nm, dm, mpopt);
        end
    end     %% methods
end         %% classdef

function TorF = ismpc2(d)
    TorF = ischar(d) || isstruct(d) && isfield(d, 'bus') && ...
        isfield(d, 'gen') && isfield(d, 'branch') && ...
        isfield(d, 'version') && strcmp(d.version, '2');
end
