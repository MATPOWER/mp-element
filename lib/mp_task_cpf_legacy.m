classdef mp_task_cpf_legacy < mp_task_cpf & tsk_shared_legacy
%MP_TASK_CPF_LEGACY  MATPOWER task for legacy continuation power flow (CPF).
%   MP_TASK_CPF_LEGACY provides implementation for continuation power flow problem.
%
%   Properties
%       ?
%
%   Methods
%       ?
%
%   See also MP_TASK_CPF

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
        function [d, mpopt] = run_pre(obj, d, mpopt)
            [d{1}, mpopt] = obj.run_pre_legacy(d{1}, mpopt);
            [d{2}, mpopt] = obj.run_pre_legacy(d{2}, mpopt);
            [d, mpopt] = run_pre@mp_task_cpf(obj, d, mpopt);
        end

        function obj = run_post(obj, mm, nm, dm, mpopt);
            if obj.nm.np ~= 0
                obj.dm.source = obj.dmc.export(obj.dm, obj.dm.source, obj.tag);
            end
        end

        %%-----  other methods  -----
        function [results, success] = legacy_post_run(obj, mpopt)
            success = obj.success;
            results = obj.dm.source;
        end
    end     %% methods
end         %% classdef
