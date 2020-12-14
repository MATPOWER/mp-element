classdef mp_task_opf < mp_task
%MP_TASK_OPF  MATPOWER task for power flow (OPF).
%   MP_TASK_OPF provides implementation for optimal power flow problem.
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

    properties
        dc      %% true if DC network model (cached, from mpopt)
    end
    
    methods
        %%-----  constructor  -----
        function obj = mp_task_opf()
            %% call parent constructor
            obj@mp_task();

            obj.tag = 'OPF';
            obj.name = 'Optimal Power Flow';
        end

        %%-----  task methods  -----
        function m = run_pre(obj, m, mpopt)
            m = run_pre@mp_task(obj, m, mpopt);     %% call parent

            %% cache some flags
            obj.dc = strcmp(upper(mpopt.model), 'DC');

            if ~obj.dc
                alg = upper(mpopt.opf.ac.solver);

                %% check for unsupported solver selection
                if strcmp(alg, 'MINOPF') || strcmp(alg, 'PDIPM') || ...
                        strcmp(alg, 'TRALM') || strcmp(alg, 'SDPOPF')
                    error('mp_task_opf/run_pre: Option ''opf.solver.ac''=''%s'' not supported.', alg);
                end
            end
        end

        %%-----  data model methods  -----
        function dm = data_model_build_post(obj, dm, mpopt)
            dm = data_model_build_post@mp_task(obj, dm, mpopt); %% call parent

            %% pre-process inputs for legacy user vars, constraints, costs
            dm.legacy_user_mod_inputs(mpopt, obj.dc);

            if ~obj.dc
                %% if requested, adjust bus voltage magnitude
                %% limits based on generator Vg setpoint
                use_vg = mpopt.opf.use_vg;
                if use_vg
                    dm.set_bus_v_lims_via_vg(use_vg);
                end
            end
        end

        function dm = data_model_update(obj, mm, nm, dm, mpopt)
            %% e.g. update data model with network model solution
            if mpopt.verbose, fprintf('-- %s data_model_update()\n', obj.tag); end
        end

        %%-----  network model methods  -----
        function nm_class = network_model_class_default(obj, dm, mpopt)
            if obj.dc
                nm_class = @mp_network_dc;
            else
                if mpopt.opf.v_cartesian
                    if mpopt.opf.current_balance
                        nm_class = @mp_network_acci;
                    else
                        nm_class = @mp_network_accs;
                    end
                else
                    if mpopt.opf.current_balance
                        nm_class = @mp_network_acpi;
                    else
                        nm_class = @mp_network_acps;
%                        nm_class = @mp_network_acps_test_nln;
                    end
                end
            end
        end

        function nm = network_model_update(obj, mm, nm)
            fprintf('-- %s network_model_update()\n', obj.tag);

            %% convert back to complex voltage vector
            x = mm.soln.x;
            out = mm.soln.output;

            if isfield(out, 'iterations')
                i = out.iterations;
            else
                i = -1;
            end
        end

        %%-----  mathematical model methods  -----
        function mm = math_model_create(obj, nm, dm, mpopt)
            %% switch back to simple opt_model, if possible
            %% I believe opf_model is required for callback functions
            %% that extract mpc from the mm
            mm = opf_model(dm.mpc);
            obj.mm = mm;
        end

        function obj = math_model_add_vars(obj, mm, nm, dm, mpopt)
            nm.add_opf_vars(nm, mm, dm, mpopt);
        end

        function obj = math_model_add_constraints(obj, mm, nm, dm, mpopt)
            nm.add_opf_constraints(nm, mm, dm, mpopt);
        end

        function obj = math_model_add_costs(obj, mm, nm, dm, mpopt)
            nm.add_opf_costs(nm, mm, dm, mpopt);
        end

        function opt = math_model_opt(obj, mm, nm, dm, mpopt)
            if strcmp(mm.problem_type(), 'NLP')
                opt = mpopt2nlpopt(mpopt, mm.problem_type(), 'DEFAULT');
            else
                opt = mpopt2qpopt(mpopt, mm.problem_type(), 'DEFAULT');
            end

            obj.mm_opt = opt;
        end

        function mm = math_model_build_post(obj, mm, nm, dm, mpopt)
            mm = dm.run_userfcn(mm, mpopt);
        end
    end     %% methods
end         %% classdef
