classdef mp_task_cpf < mp_task_pf
%MP_TASK_CPF  MATPOWER task for continuation power flow (CPF).
%   MP_TASK_CPF provides implementation for continuation power flow problem.
%
%   Properties
%       ?
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

%     properties
%     end
    
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

        %%-----  data model methods  -----

        %%-----  network model methods  -----

        %%-----  mathematical model methods  -----
        function mm = math_model_build_pre(obj, mm, nm, dm, mpopt)
            mm.userdata.aux_data = nm.cpf_aux_data(dm, mpopt);
        end

        function obj = math_model_add_vars(obj, mm, nm, dm, mpopt)
            nm.cpf_add_vars(mm, nm, dm, mpopt);
        end

        function obj = math_model_add_constraints(obj, mm, nm, dm, mpopt)
            nm.cpf_add_constraints(mm, nm, dm, mpopt);
        end

        function opt = math_model_opt(obj, mm, nm, dm, mpopt)
            opt = nm.cpf_solve_opts(mm, dm, mpopt);
            obj.mm_opt = opt;
        end
    end     %% methods
end         %% classdef
