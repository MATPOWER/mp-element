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

        function dm = next_dm(obj, mm, nm, dm, mpopt)
            if ~obj.dc && mpopt.pf.enforce_q_lims
                %% adjust iteration count for previous runs
                obj.iterations = obj.iterations + mm.soln.output.iterations;
                mm.soln.output.iterations = obj.iterations;

                %% enforce Q limits
                [success, dm] = obj.enforce_q_lims(nm, dm, mpopt);
                if ~success                 %% entire task fails if Q lim
                    obj.success = success;  %% enforcement indicates failure
                end
            else        %% don't enforce generator Q limits, once is enough
                dm = [];
            end
        end

        function [success, dm] = enforce_q_lims(obj, nm, dm, mpopt);
            gen_dme = dm.elements.gen;
            [mn, mx, both] = gen_dme.violated_q_lims(dm, mpopt);

            if ~isempty(both)   %% we have some Q limit violations
                if isempty(mn) && isempty(mx)   %% infeasible
                    if mpopt.verbose
                        fprintf('All %d remaining gens exceed their Q limits : INFEASIBLE PROBLEM\n', length(both));
                    end
                    dm = [];
                    success = 0;
                else
                    if mpopt.verbose && ~isempty(mx)
                        fprintf('Gen %d at upper Q limit, converting to PQ bus\n', gen_dme.on(mx));
                    end
                    if mpopt.verbose && ~isempty(mn)
                        fprintf('Gen %d at lower Q limit, converting to PQ bus\n', gen_dme.on(mn));
                    end

                    %% save corresponding limit values
                    obj.fixed_q_qty(mx) = gen_dme.qg_ub(mx);
                    obj.fixed_q_qty(mn) = gen_dme.qg_lb(mn);
                    mx = [mx;mn];

                    %% set qg to binding limit
                    gen_dme.tab.qg(gen_dme.on(mx)) = ...
                        obj.fixed_q_qty(mx) * dm.base_mva;

                    %% convert to PQ bus
                    bus_dme = dm.elements.bus;
                    ref0 = find(bus_dme.type == NODE_TYPE.REF);
                    bidx = bus_dme.i2on(gen_dme.bus(gen_dme.on(mx)));   %% bus of mx
                    if length(ref0) > 1 && any(bus_dme.type(bidx) == NODE_TYPE.REF)
                        error('mp_data/enforce_q_lims: Sorry, MATPOWER cannot enforce Q limits for slack buses in systems with multiple slacks.');
                    end
                    %% set bus type to PQ
                    bus_dme.set_bus_type_pq(dm, bidx);
                    %% potentially pick new slack bus
                    ntv = nm.node_types(nm, dm);        %% node type vector
                    [i1, iN] = nm.get_node_idx('bus');  %% bus node indices
                    btv = ntv(i1:iN);                   %% bus type vector

                    %% indicate if there's been a change in slack bus
                    ref = find(btv == NODE_TYPE.REF);   %% new ref bus indices
                    if mpopt.verbose && ref ~= ref0
                        fprintf('Bus %d is new slack bus\n', ...
                            bus_dme.ID(bus_dme.on(ref)));
                    end

                    %% save indices to list of Q limited gens
                    obj.fixed_q_idx = [obj.fixed_q_idx; mx];

                    %% update dm for next step
                    dm.initialize();
                    dm.update_status();
                    dm.build_params();
                    success = 1;
                end
            else                %% no more Q violations
                dm = [];
                success = 1;
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
        function mm_class = math_model_class_default(obj, nm, dm, mpopt)
            switch upper(mpopt.model)
                case 'AC'
                    if mpopt.pf.v_cartesian
                        if mpopt.pf.current_balance
                            mm_class = @mp_math_pf_acci;
                        else
                            mm_class = @mp_math_pf_accs;
                        end
                    else
                        if mpopt.pf.current_balance
                            mm_class = @mp_math_pf_acpi;
                        else
                            mm_class = @mp_math_pf_acps;
                        end
                    end
                case 'DC'
                    mm_class = @mp_math_pf_dc;
            end
        end
    end     %% methods
end         %% classdef
