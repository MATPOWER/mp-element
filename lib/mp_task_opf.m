classdef mp_task_opf < mp_task
%MP_TASK_OPF  MATPOWER task for optimal power flow (OPF).
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
        function dm = data_model_build_post(obj, dm, dmc, mpopt)
            %% call parent
            dm = data_model_build_post@mp_task(obj, dm, dmc, mpopt);

            if ~obj.dc
                %% if requested, adjust bus voltage magnitude
                %% limits based on generator vm_setpoint
                use_vg = mpopt.opf.use_vg;
                if use_vg
                    dm.set_bus_v_lims_via_vg(use_vg);
                end
            end
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

        %%-----  mathematical model methods  -----
        function mm_class = math_model_class(obj, nm, dm, mpopt)
            mm_class = @mp_math_opf;
        end
    end     %% methods
end         %% classdef
