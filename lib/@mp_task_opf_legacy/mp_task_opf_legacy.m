classdef mp_task_opf_legacy < mp_task_opf & tsk_shared_legacy
%MP_TASK_OPF_LEGACY  MATPOWER task for legacy optimal power flow (OPF).
%   MP_TASK_OPF_LEGACY provides implementation for optimal power flow problem.
%
%   Properties
%       ?
%
%   Methods
%       ?
%
%   See also MP_TASK_OPF

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
    end

    methods
        %%-----  task methods  -----
        function [d, mpopt] = run_pre(obj, d, mpopt)
            [d, mpopt] = obj.run_pre_legacy(d, mpopt);
            [d, mpopt] = run_pre@mp_task_opf(obj, d, mpopt);
        end

        function obj = run_post(obj, mm, nm, dm, mpopt);
            if obj.nm.np ~= 0
                obj.dm.source = obj.dmc.export(obj.dm, obj.dm.source);
            end
        end

        %%-----  data model converter methods  -----
        function dmc_class = dm_converter_class_mpc2_default(obj)
            dmc_class = @mp_dm_converter_mpc2_legacy;
        end

        %%-----  data model methods  -----
        function dm = data_model_build_post(obj, dm, dmc, mpopt)
            %% call parent
            dm = data_model_build_post@mp_task_opf(obj, dm, dmc, mpopt);

            %% pre-process inputs for legacy user vars, constraints, costs
            dm = dmc.legacy_user_mod_inputs(dm, mpopt, obj.dc);
        end

        %%-----  mathematical model methods  -----
        function mm_class = math_model_class_default(obj, nm, dm, mpopt)
            %% mm_shared_opf_legacy (compatible with opf_model) is required
            %% to support legacy cost functions and callback functions that
            %% expect to find mpc in mm.mpc.
            switch upper(mpopt.model)
                case 'AC'
                    if mpopt.opf.v_cartesian
                        if mpopt.opf.current_balance
                            mm_class = @mp_math_opf_acci_legacy;
                        else
                            mm_class = @mp_math_opf_accs_legacy;
                        end
                    else
                        if mpopt.opf.current_balance
                            mm_class = @mp_math_opf_acpi_legacy;
                        else
                            mm_class = @mp_math_opf_acps_legacy;
                        end
                    end
                case 'DC'
                    mm_class = @mp_math_opf_dc_legacy;
            end
        end
    end     %% methods
end         %% classdef
