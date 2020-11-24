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
%
%   Methods
%       setup() - sets up the mathematical model
%       network_model() - updates and returns the network model with the
%           current mathematical model solution

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
    end

    methods
        function dm_class = data_model_class(obj, m, mpopt)
            dm_class = @mp_data_mpc2;
        end

        function dm = create_data_model(obj, m, mpopt)
            if ~isa(m, 'mp_data')
                dm_class = obj.data_model_class(m, mpopt);
                dm = dm_class(m);
            else
                dm = m;
            end
            obj.dm = dm;
        end

%         nm_class = network_model_class(obj, dm, mpopt)

        function nm = create_network_model(obj, dm, mpopt)
            nm_class = obj.network_model_class(dm, mpopt);
            nm = nm_class();
            nm.create_model(dm, mpopt);
            obj.nm = nm;
        end

        function mm_class = math_model_class(obj, nm, dm, mpopt)
            mm_class = @opt_model;
        end

        function mm = create_math_model(obj, nm, dm, mpopt)
            mm_class = obj.math_model_class(nm, dm, mpopt);
            mm = mm_class();

            %% construct auxiliary data
            obj.add_aux_data(mm, nm, dm, mpopt);

            %% add variables, constraints, costs
            if nm.np ~= 0       %% skip for empty model
                obj.add_vars(mm, nm, dm, mpopt);
                obj.add_constraints(mm, nm, dm, mpopt);
                obj.add_costs(mm, nm, dm, mpopt);
                
                %% add user customization
                mm = obj.customize_mm(mm, nm, dm, mpopt);
            end
            obj.mm = mm;
        end

        function obj = add_aux_data(obj, mm, nm, dm, mpopt)
        end

        function obj = add_vars(obj, mm, nm, dm, mpopt)
        end

        function obj = add_constraints(obj, mm, nm, dm, mpopt)
        end

        function obj = add_costs(obj, mm, nm, dm, mpopt)
        end

        function mm = customize_mm(obj, mm, nm, dm, mpopt)
            %% execute userfcn callbacks for 'formulation' stage
            mpc = dm.mpc;
            if isfield(mpc, 'userfcn')
                userfcn = mpc.userfcn;
            else
                userfcn = [];
            end
            mm = run_userfcn(userfcn, 'formulation', mm, mpopt);
        end

        function obj = add_mm_opt(obj, mpopt)
        end

        function success = run(obj, mm, nm, dm, mpopt)
            %% get solve options
            obj.add_mm_opt(mpopt);

            %% solve mathematical model
            if obj.mm_opt.verbose, fprintf('-----  mp_task_pf.run()  -----\n'); end
            mm.solve(obj.mm_opt);
            success = (mm.soln.eflag > 0);
        end

        function print_soln(obj, fname)
            if isempty(fname)
                fprintf('-- mp_task/print_soln()\n');
            else
                fprintf('-- mp_task/print_soln(''%s'')\n', fname);
            end
        end

        function save_soln(obj, fname)
            fprintf('-- mp_task/save_soln(''%s'')\n', fname);
        end
    end     %% methods
end         %% classdef
