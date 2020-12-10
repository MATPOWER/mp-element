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
%       success - success flag, 1 - math model solved, 0 - didn't solve
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
        dm      %% data model object
        nm      %% network model object
        mm      %% mathematical model object
        mm_opt  %% solve options for mathematical model
        tag     %% task tag - e.g. 'PF', 'CPF', 'OPF'
        name    %% task name - e.g. 'Power Flow', etc.
        success %% success flag, 1 - math model solved, 0 - didn't solve
    end

    methods
        function obj = run(obj, m, mpopt)
            %% create models
            dm = obj.data_model_build(m, mpopt);
            nm = obj.network_model_build(dm, mpopt);
            mm = obj.math_model_build(nm, dm, mpopt);

            %% get solve options
            mm_opt = obj.math_model_opt(mm, nm, dm, mpopt);

            %% solve mathematical model
            if obj.mm_opt.verbose
                fprintf('-----  SOLVE %s  -----\n', obj.tag);
            end
            mm.solve(mm_opt);
            obj.success = (mm.soln.eflag > 0);

%             if obj.success
                %% update network model with math model solution
                nm = obj.network_model_update(mm, nm);

                %% update data model with network model solution
                dm = obj.data_model_update(mm, nm, dm, mpopt);
                dm = obj.data_model_update_post(mm, nm, dm, mpopt);
%             end
        end

        function print_soln(obj, fname)
            if isempty(fname)
                fprintf('-- %s print_soln()\n', obj.tag);
            else
                fprintf('-- %s print_soln(''%s'')\n', obj.tag, fname);
            end
        end

        function save_soln(obj, fname)
            fprintf('-- %s save_soln(''%s'')\n', obj.tag, fname);
        end

        %%-----  data model methods  -----
        function dm_class = data_model_class(obj, m, mpopt)
            if isfield(mpopt.exp, 'data_model_class') && ...
                    ~isempty(mpopt.exp.data_model_class)
                dm_class = mpopt.exp.data_model_class;
            else
                dm_class = @mp_data_mpc2;
            end
        end

        function dm = data_model_build(obj, m, mpopt)
            if ~isa(m, 'mp_data')
                dm_class = obj.data_model_class(m, mpopt);
                if isfield(mpopt.exp, 'dm_element_classes')
                    dm = dm_class(m, mpopt.exp.dm_element_classes);
                else
                    dm = dm_class(m);
                end
            else
                dm = m;
            end
            obj.dm = dm;
            dm = obj.data_model_build_post(dm, mpopt);
        end

        function dm = data_model_build_post(obj, dm, mpopt)
%             dm.ext2int(mpopt);
        end

%         function dm = data_model_update(obj, mm, nm, dm, mpopt)
%             %% e.g. update data model with network model solution
%         end

        function dm = data_model_update_post(obj, mm, nm, dm, mpopt)
%             if mpopt.verbose, fprintf('-- %s data_model_update_post()\n', obj.tag); end
%             dm.int2ext(mpopt);
        end

        %%-----  network model methods  -----
        function nm_class = network_model_class_override(obj, dm, mpopt)
            if isfield(mpopt.exp, 'network_model_class') && ...
                    ~isempty(mpopt.exp.network_model_class)
                nm_class = mpopt.exp.network_model_class;
            else
                nm_class = [];
            end
        end

        function nm = network_model_build(obj, dm, mpopt)
            nm_class = obj.network_model_class(dm, mpopt);
            nm = nm_class();
            obj.nm = nm;

            nm = obj.network_model_build_pre(nm, dm, mpopt);
            nm.build(dm, mpopt);
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

        function nm = network_model_update(obj, mm, nm)
        end

        %%-----  mathematical model methods  -----
        function mm_class = math_model_class(obj, nm, dm, mpopt)
            mm_class = @opt_model;
        end

        function mm = math_model_build(obj, nm, dm, mpopt)
            mm_class = obj.math_model_class(nm, dm, mpopt);
            mm = mm_class();
            obj.mm = mm;

            mm = obj.math_model_build_pre(mm, nm, dm, mpopt);

            %% add variables, constraints, costs
            if nm.np ~= 0       %% skip for empty model
                obj.math_model_add_vars(mm, nm, dm, mpopt);
                obj.math_model_add_constraints(mm, nm, dm, mpopt);
                obj.math_model_add_costs(mm, nm, dm, mpopt);
                
                %% add user customization
                mm = obj.math_model_build_post(mm, nm, dm, mpopt);
            end
        end

        function mm = math_model_build_pre(obj, mm, nm, dm, mpopt)
        end

        function mm = math_model_build_post(obj, mm, nm, dm, mpopt)
        end

        function obj = math_model_add_vars(obj, mm, nm, dm, mpopt)
        end

        function obj = math_model_add_constraints(obj, mm, nm, dm, mpopt)
        end

        function obj = math_model_add_costs(obj, mm, nm, dm, mpopt)
        end

        function opt = math_model_opt(obj, mm, nm, dm, mpopt)
        end
    end     %% methods
end         %% classdef
