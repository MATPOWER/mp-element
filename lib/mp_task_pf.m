classdef mp_task_pf < mp_task
%MP_TASK_PF  MATPOWER task for power flow (PF).
%   MP_TASK_PF provides implementation for power flow problem.
%
%   Properties
%       dc
%       iterations
%       ref
%       ref0
%       va_ref0
%       fixed_q_idx
%       fixed_q_qty
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2020-2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        dc      %% true if DC network model (cached in run_pre(), from mpopt)
        iterations
        ref
        ref0
        va_ref0
        fixed_q_idx
        fixed_q_qty
    end

    methods
        %%-----  constructor  -----
        function obj = mp_task_pf()
            %% call parent constructor
            obj@mp_task();

            obj.tag = 'PF';
            obj.name = 'Power Flow';
        end

        %%-----  task methods  -----
        function [d, mpopt] = run_pre(obj, d, mpopt)
            [d, mpopt] = run_pre@mp_task(obj, d, mpopt);    %% call parent

            %% cache DC model flag
            obj.dc = strcmp(upper(mpopt.model), 'DC');
        end

        function [results, success] = legacy_post_run(obj, mpopt)
            success = obj.success;
            results = obj.dm.mpc;
            if obj.nm.np ~= 0
                results = obj.dmc.elements.bus.export( ...
                    obj.dm.elements.bus, results, ...
                    {'type', 'vm', 'va'});
                results = obj.dmc.elements.gen.export( ...
                    obj.dm.elements.gen, results, ...
                    {'pg', 'qg'});
                results = obj.dmc.elements.branch.export( ...
                    obj.dm.elements.branch, results, ...
                    {'ql_to', 'pl_to', 'ql_fr', 'pl_fr'});
%                 results = obj.dmc.export(obj.dm, obj.dm.mpc);
            end
        end

        function dm = next_dm(obj, mm, nm, dm, mpopt)
            if ~obj.dc && mpopt.pf.enforce_q_lims
                %% adjust iteration count for previous runs
                obj.iterations = obj.iterations + mm.soln.output.iterations;
                mm.soln.output.iterations = obj.iterations;

                %% enforce Q limits
                [success, d, obj] = dm.pf_enforce_q_lims(obj, nm, mpopt);
                if ~success                 %% entire task fails if Q lim
                    obj.success = success;  %% enforcement indicates failure
                end
                if isempty(d)   %% Q limits are satisfied (or failed)
                    dm = [];
                else            %% use new data model to satisfy limits
                    dm = obj.data_model_build(d, obj.dmc, mpopt);
                end
            else        %% don't enforce generator Q limits, once is enough
                dm = [];
            end
        end

        %%-----  data model methods  -----

        %%-----  network model methods  -----
        function nm_class = network_model_class_default(obj, dm, mpopt)
            switch upper(mpopt.model)
                case 'AC'
                    if mpopt.pf.v_cartesian
                        if mpopt.pf.current_balance
                            nm_class = @mp_network_acci;
                        else
                            nm_class = @mp_network_accs;
%                                nm_class = @mp_network_accs_test_nln;
                        end
                    else
                        if mpopt.pf.current_balance
                            nm_class = @mp_network_acpi;
                        else
                            nm_class = @mp_network_acps;
%                                nm_class = @mp_network_acps_test_nln;
                        end
                    end
                case 'DC'
                    nm_class = @mp_network_dc;
            end
        end

        function nm = network_model_build_post(obj, nm, dm, mpopt)
            %% initialize task data, if non-empty AC case with Q lim enforced
            if ~obj.dc && mpopt.pf.enforce_q_lims ~= 0 && nm.np ~= 0
                if obj.i_nm == 1
                    [ref, ~, ~] = nm.node_types(obj, dm);
                    gen_dme =  dm.elements.gen;
                    obj.iterations = 0;
                    obj.ref0 = ref;             %% initial ref node indices
                    obj.ref = ref;              %% current ref node indices
                    obj.va_ref0 = nm.get_va(ref);%% initial ref node voltage angles
                    obj.fixed_q_idx = [];       %% indices of fixed Q gens
                    obj.fixed_q_qty = zeros(gen_dme.n, 1);  %% Q output of fixed Q gens
                else        %% update index of ref bus
                    [obj.ref, ~, ~] = nm.node_types(obj, dm);
                end
            end
        end

        function nm = network_model_x_soln(obj, mm, nm)
            nm = network_model_x_soln@mp_task(obj, mm, nm);

            %% if ref node has been changed, adjust voltage angles
            %% to make angle at original ref node = specified value
            if ~obj.dc && obj.i_nm > 1 && obj.ref ~= obj.ref0
                vm = abs(nm.soln.v);
                va = angle(nm.soln.v);
                va = va - va(obj.ref0) + obj.va_ref0;
                nm.soln.v = vm .* exp(1j * va);
            end
        end

        %%-----  mathematical model methods  -----
        function mm_class = math_model_class(obj, nm, dm, mpopt)
            mm_class = @mp_math_pf;
        end

        function opt = math_model_opt(obj, mm, nm, dm, mpopt)
            opt = math_model_opt@mp_task(obj, mm, nm, dm, mpopt);

            %%-----  HACK ALERT  -----
            %% only needed for fast-decoupled power flow,
            %% used by mp_data_mpc2/fdpf_B_matrix_models()
            dm.userdata.dmc = obj.dmc;
            %%-----  end of HACK  -----
        end
    end     %% methods
end         %% classdef
