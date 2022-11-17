classdef (Abstract) task < handle
%MP.TASK  MATPOWER task abstract base class.
%   Each task type (e.g. power flow, CPF, OPF) will inherit from MP.TASK.
%
%   MP.TASK provides properties and methods related to the specific
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
%       et - elapsed time (seconds) for run() method
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

    properties (Abstract)
        tag     %% task tag - e.g. 'PF', 'CPF', 'OPF'
        name    %% task name - e.g. 'Power Flow', etc.
    end
    properties
        dmc     %% data model converter object
        dm      %% data model object
        nm      %% network model object
        mm      %% mathematical model object
        mm_opt  %% solve options for mathematical model
        i_dm    %% iteration counter for data model loop
        i_nm    %% iteration counter for network model loop
        i_mm    %% iteration counter for math model loop
        success %% success flag, 1 - math model solved, 0 - didn't solve
        message %% output message
        et      %% elapsed time (seconds) for run() method
    end

    methods
        %%-----  task methods  -----
        function obj = run(obj, d, mpopt, mpx)
            t0 = tic;       %% start timer
            if nargin < 4
                mpx = {};   %% no MATPOWER extensions by default
            end

            [d, mpopt] = obj.run_pre(d, mpopt);

            dmc = obj.dm_converter_create(d, mpopt, mpx);
            obj.dmc = dmc;

            %% initialize
            obj.i_dm = 0;   %% iteration counter for data model loop
            obj.i_nm = 0;   %% iteration counter for network model loop
            obj.i_mm = 0;   %% iteration counter for math model loop

            %% build data model
            obj.i_dm = obj.i_dm + 1;
            dm = obj.data_model_build(d, dmc, mpopt, mpx);
            obj.dm = dm;    %% stash current data model in task object

            while ~isempty(dm)  %% begin data model loop
                %% build network model
                obj.i_nm = obj.i_nm + 1;
                nm = obj.network_model_build(dm, mpopt, mpx);
                obj.nm = nm;    %% stash current network model in task object

                while ~isempty(nm)  %% begin network model loop
                    %% build math model
                    obj.i_mm = obj.i_mm + 1;
                    mm = obj.math_model_build(nm, dm, mpopt, mpx);
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

                        [mm, nm, dm] = obj.next_mm(mm, nm, dm, mpopt, mpx);
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

                        [nm, dm] = obj.next_nm(mm, nm, dm, mpopt, mpx);
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

                dm = obj.next_dm(mm, nm, dm, mpopt, mpx);
                if ~isempty(dm)
                    obj.dm = dm;
                    obj.i_dm = obj.i_dm + 1;
                end
            end                 %% end data model loop
            obj.run_post(mm, nm, obj.dm, mpopt);
            obj.et = toc(t0);   %% stop timer
        end

        function [mm, nm, dm] = next_mm(obj, mm, nm, dm, mpopt, mpx)
            %% return new math model, or empty matrix if finished
            mm = [];
        end

        function [nm, dm] = next_nm(obj, mm, nm, dm, mpopt, mpx)
            %% return new network model, or empty matrix if finished
            nm = [];
        end

        function dm = next_dm(obj, mm, nm, dm, mpopt, mpx)
            %% return new data model, or empty matrix if finished
            dm = [];
        end

        function [d, mpopt] = run_pre(obj, d, mpopt)
        end

        function obj = run_post(obj, mm, nm, dm, mpopt);
        end

        function print_soln(obj, mpopt, fname)
            if nargin < 3
                fname = '';
            end

            %% print to a file, if requested
            if fname
                [fd, msg] = fopen(fname, 'at');
                if fd == -1
                    warning('mp.task/print_soln: could not open file ''%s'' for writing\n%s', fname, msg);
                else
                    obj.print_soln_header(mpopt, fd);
                    if mpopt.out.all == 0
                        obj.dm.pretty_print(mpoption(mpopt, 'out.all', -1), fd);
                    else
                        obj.dm.pretty_print(mpopt, fd);
                    end
                    fclose(fd);
                end
            end

            %% print to stdio
            if mpopt.out.all
                obj.print_soln_header(mpopt);
                if obj.success || mpopt.out.force
                    obj.dm.pretty_print(mpopt);
                end
            end
        end

        function print_soln_header(obj, mpopt, fd)
            if nargin < 3
                fd = 1;     %% print to stdio by default
            end

            %% succeeded/failed line
            if obj.success
                succ_fail = 'succeeded';
            else
                succ_fail = 'failed';
            end
            fprintf(fd, ...
                '\n%s %s in %.2f seconds (%.2f setup + %.2f solve)\n', ...
                obj.tag, succ_fail, obj.et, obj.et - obj.mm.soln.output.et, ...
                obj.mm.soln.output.et);
        end

        function save_soln(obj, fname)
            %% export solution
            if obj.nm.np ~= 0
                obj.dm.source = obj.dmc.export(obj.dm, obj.dm.source);
            end

            %% save exported solution
            obj.dmc.save(fname, obj.dm.source);
        end

        %%-----  data model converter methods  -----
        function dmc_class = dm_converter_class(obj, d, mpopt, mpx)
            %% manual override
            if isfield(mpopt.exp, 'dm_converter_class') && ...
                    ~isempty(mpopt.exp.dm_converter_class)
                dmc_class = mpopt.exp.dm_converter_class;
            else    %% use default
                %% detect input data type
                if ismpc2(d)
                    d_fmt = 'mpc2';
                else
                    error('mp.task: input data format not recognized');
                end

                %% get default class
                switch d_fmt
                    case 'mpc2'
                        dmc_class = obj.dm_converter_class_mpc2_default();
                    otherwise
                        error('mp.task: input data format not recognized');
                end

                %% apply extensions
                for k = 1:length(mpx)
                    dmc_class = mpx{k}.dm_converter_class(dmc_class, d_fmt, mpopt);
                end
            end
        end

        function dmc_class = dm_converter_class_mpc2_default(obj)
            dmc_class = @mp.dm_converter_mpc2;
        end

        function dmc = dm_converter_create(obj, d, mpopt, mpx)
            if isa(d, 'mp.data_model')
                dmc = [];
            else
                dmc_class = obj.dm_converter_class(d, mpopt, mpx);
                dmc = dmc_class();

                %% apply extensions
                for k = 1:length(mpx)
                    dmc_elements = mpx{k}.dmc_element_classes(dmc_class, ...
                                                    dmc.format_tag, mpopt);
                    if ~isempty(dmc_elements)
                        dmc.modify_element_classes(dmc_elements);
                    end
                end

                %% apply user-supplied dmc.element_classes overrides
                if isfield(mpopt.exp, 'dmc_element_classes') && ...
                        ~isempty(mpopt.exp.dmc_element_classes)
                    dmc.modify_element_classes(mpopt.exp.dmc_element_classes);
                end

                dmc.build();

                %% remove excluded elements (results in corresponding elements
                %% being excluded from data, network and math models as well)
                if isfield(mpopt.exp, 'exclude_elements') && ...
                        ~isempty(mpopt.exp.exclude_elements)
                    ex = mpopt.exp.exclude_elements;
                    for k = length(ex):-1:1
                        if ~dmc.elements.is_index_name(ex{k})
                            ex(k) = [];     %% skip missing exclusions
                        end
                    end
                    dmc.elements.delete_elements(ex);
                end
            end
        end

        %%-----  data model methods  -----
        function dm_class = data_model_class(obj, d, mpopt, mpx)
            %% manual override
            if isfield(mpopt.exp, 'data_model_class') && ...
                    ~isempty(mpopt.exp.data_model_class)
                dm_class = mpopt.exp.data_model_class;
            else    %% use default
                %% get default class
                dm_class = obj.data_model_class_default();

                %% apply extensions
                for k = 1:length(mpx)
                    dm_class = mpx{k}.data_model_class(dm_class, obj.tag, mpopt);
                end
            end
        end

        function dm_class = data_model_class_default(obj)
            dm_class = @mp.data_model;
        end

        function dm = data_model_create(obj, d, mpopt, mpx)
            dm_class = obj.data_model_class(d, mpopt, mpx);
            dm = dm_class();

            %% apply extensions
            for k = 1:length(mpx)
                dm_elements = mpx{k}.dm_element_classes(dm_class, obj.tag, mpopt);
                if ~isempty(dm_elements)
                    dm.modify_element_classes(dm_elements);
                end
            end

            %% apply user-supplied dm.element_classes overrides
            if isfield(mpopt.exp, 'dm_element_classes') && ...
                    ~isempty(mpopt.exp.dm_element_classes)
                dm.modify_element_classes(mpopt.exp.dm_element_classes);
            end
        end

        function dm = data_model_build(obj, d, dmc, mpopt, mpx)
            if isa(d, 'mp.data_model')
                dm = d;
            else
                dm = obj.data_model_create(d, mpopt, mpx);
                [dm, d] = obj.data_model_build_pre(dm, d, dmc, mpopt);
                dm.build(d, dmc);
                dm = obj.data_model_build_post(dm, dmc, mpopt);
            end
        end

        function [dm, d] = data_model_build_pre(obj, dm, d, dmc, mpopt)
        end

        function dm = data_model_build_post(obj, dm, dmc, mpopt)
        end

        %%-----  network model methods  -----
        function nm_class = network_model_class(obj, dm, mpopt, mpx)
            %% manual override
            if isfield(mpopt.exp, 'network_model_class') && ...
                    ~isempty(mpopt.exp.network_model_class)
                nm_class = mpopt.exp.network_model_class;
            else    %% use default
                %% get default class
                nm_class = obj.network_model_class_default(dm, mpopt);

                %% apply extensions
                for k = 1:length(mpx)
                    nm_class = mpx{k}.network_model_class(nm_class, obj.tag, mpopt);
                end
            end
        end

        function nm_class = network_model_class_default(obj, dm, mpopt)
            error('mp.task/network_model_class_default: must be implemented in subclass');
        end

        function nm = network_model_create(obj, dm, mpopt, mpx)
            nm_class = obj.network_model_class(dm, mpopt, mpx);
            nm = nm_class();
            nm.init_set_types();

            %% apply extensions
            for k = 1:length(mpx)
                nm_elements = mpx{k}.nm_element_classes(nm_class, obj.tag, mpopt);
                if ~isempty(nm_elements)
                    nm.modify_element_classes(nm_elements);
                end
            end

            %% apply user-supplied nm.element_classes overrides
            if isfield(mpopt.exp, 'nm_element_classes') && ...
                    ~isempty(mpopt.exp.nm_element_classes)
                nm.modify_element_classes(mpopt.exp.nm_element_classes);
            end
        end

        function nm = network_model_build(obj, dm, mpopt, mpx)
            nm = obj.network_model_create(dm, mpopt, mpx);
            nm = obj.network_model_build_pre(nm, dm, mpopt);
            nm.build(dm);
            nm = obj.network_model_build_post(nm, dm, mpopt);
        end

        function nm = network_model_build_pre(obj, nm, dm, mpopt)
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
        function mm_class = math_model_class(obj, nm, dm, mpopt, mpx)
            %% manual override
            if isfield(mpopt.exp, 'math_model_class') && ...
                    ~isempty(mpopt.exp.math_model_class)
                mm_class = mpopt.exp.math_model_class;
            else    %% use default
                %% get default class
                mm_class = obj.math_model_class_default(nm, dm, mpopt);

                %% apply extensions
                for k = 1:length(mpx)
                    mm_class = mpx{k}.math_model_class(mm_class, obj.tag, mpopt);
                end
            end
        end

        function mm_class = math_model_class_default(obj, nm, dm, mpopt)
            error('mp.task/math_model_class_default: must be implemented in subclass');
        end

        function mm = math_model_create(obj, nm, dm, mpopt, mpx)
            mm_class = obj.math_model_class(nm, dm, mpopt, mpx);
            mm = mm_class();
            mm.init_set_types();

            %% apply extensions
            for k = 1:length(mpx)
                mm_elements = mpx{k}.mm_element_classes(mm_class, obj.tag, mpopt);
                if ~isempty(mm_elements)
                    mm.modify_element_classes(mm_elements);
                end
            end

            %% apply user-supplied mm.element_classes overrides
            if isfield(mpopt.exp, 'mm_element_classes') && ...
                    ~isempty(mpopt.exp.mm_element_classes)
                mm.modify_element_classes(mpopt.exp.mm_element_classes);
            end
        end

        function mm = math_model_build(obj, nm, dm, mpopt, mpx)
            mm = obj.math_model_create(nm, dm, mpopt, mpx);

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
