classdef mp_task_opf_legacy < mp_task_opf
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
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
    end

    methods
        %%-----  task methods  -----
        function obj = run_post(obj, mm, nm, dm, mpopt);
            if obj.nm.np ~= 0
                obj.dm.userdata.mpc = ...
                    obj.dmc.export(obj.dm, obj.dm.userdata.mpc, obj.tag);
            end
        end

        %%-----  data model methods  -----
        function dm = data_model_build_post(obj, dm, dmc, mpopt)
            %% call parent
            dm = data_model_build_post@mp_task_opf(obj, dm, dmc, mpopt);

            %% pre-process inputs for legacy user vars, constraints, costs
            dm = dmc.legacy_user_mod_inputs(dm, mpopt, obj.dc);
        end

        %%-----  mathematical model methods  -----
        function mm_class = math_model_class(obj, nm, dm, mpopt)
            %% mp_math_opf_legacy (based on opf_model) is required to
            %% support legacy cost functions and callback functions that
            %% expect to find mpc in mm.mpc.
            mm_class = @mp_math_opf_legacy;
        end
    end     %% methods
end         %% classdef
