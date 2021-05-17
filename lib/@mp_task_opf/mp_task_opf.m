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
        dc      %% true if DC network model (cached in run_pre(), from mpopt)
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
        function [d, mpopt] = run_pre(obj, d, mpopt)
            [d, mpopt] = run_pre@mp_task(obj, d, mpopt);     %% call parent

            %% cache DC model flag
            obj.dc = strcmp(upper(mpopt.model), 'DC');

            %% check for unsupported AC OPF solver selection
            if ~obj.dc
                alg = upper(mpopt.opf.ac.solver);
                switch alg
                    case 'IPOPT'
                        if ~have_feature('ipopt')
                            error('mp_task_opf/run_pre: MPOPT.opf.ac.solver = ''%s'' requires IPOPT (see https://github.com/coin-or/Ipopt)', alg);
                        end
                    case 'FMINCON'
                        if ~have_feature('fmincon')
                            error('mp_task_opf/run_pre: MPOPT.opf.ac.solver = ''%s'' requires FMINCON (Optimization Toolbox 2.x or later)', alg);
                        end
                    case 'KNITRO'
                        if ~have_feature('knitro')
                            error('mp_task_opf/run_pre: MPOPT.opf.ac.solver = ''%s'' requires Artelys Knitro (see https://www.artelys.com/solvers/knitro/)', alg);
                        end
                    case {'MINOPF', 'PDIPM', 'TRALM', 'SDPOPF'}
                        error('mp_task_opf/run_pre: MPOPT.opf.ac.solver = ''%s'' not supported.', alg);
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
            nm.opf_data_model_update(mm, nm, dm, mpopt);
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

        function nm = network_model_convert_x(obj, mm, nm)
            %% convert solved state from math model to network model soln
            [nm.soln.v, nm.soln.z, nm.soln.x] = nm.opf_convert_x(mm.soln.x, ...
                                                mm.get_userdata('aux_data'));
        end

        %%-----  mathematical model methods  -----
        function mm = math_model_create(obj, nm, dm, mpopt)
            %% switch back to simple opt_model, if possible
            %% I believe opf_model is required for callback functions
            %% that extract mpc from the mm
            mm = opf_model(dm.mpc).init_set_types();
            obj.mm = mm;
        end

        function obj = math_model_build_it(obj, mm, nm, dm, mpopt)
            nm.opf_build_math_model(mm, dm, mpopt);
        end

        function opt = math_model_opt(obj, mm, nm, dm, mpopt)
            opt = nm.opf_solve_opts(mm, dm, mpopt);
            obj.mm_opt = opt;
        end

        function mm = math_model_build_post(obj, mm, nm, dm, mpopt)
            if obj.dc
                %% user data required by toggle_softlims
                branch_nme = nm.elm_by_name('branch');
                [Bbr, pbr] = branch_nme.get_params(1:branch_nme.nk, {'B', 'p'});
                mm.userdata.Bf = Bbr * branch_nme.C';
                mm.userdata.Pfinj = pbr;
            end

            mm = dm.run_userfcn(mm, mpopt);
        end
    end     %% methods
end         %% classdef
