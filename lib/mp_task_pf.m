classdef mp_task_pf < mp_task
%MP_TASK_PF  MATPOWER task for power flow (PF).
%   MP_TASK_PF provides implementation for power flow problem.
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

        function [dm, td] = next_dm(obj, mm, nm, dm, td, mpopt)
            if ~obj.dc && mpopt.pf.enforce_q_lims
                %% initialize task data
                if isempty(td)
                    bus_dme =  dm.elm_by_name('bus');
                    gen_dme =  dm.elm_by_name('gen');
                    ref = find(bus_dme.isref);
                    td = struct( ...        %% task data
                        'ref', ref, ...                 %% ref bus indices
                        'Varef', bus_dme.Va0(ref), ...  %% ref bus V angles
                        'limited', [], ...              %% indices of fixed Q gens
                        'fixedQg', zeros(gen_dme.n, 1));%% Q output of fixed Q gens
                end
                [success, d, td] = dm.pf_enforce_q_lims(td, nm, mpopt);
                if ~success                 %% entire task fails if Q lim
                    obj.success = success;  %% enforcement indicates failure
                end
                if isempty(d)   %% Q limits are satisfied
                    dm = [];
                else            %% use new data model to satisfy limits
                    dm = obj.data_model_build(d, mpopt);
                end
            else        %% don't enforce generator Q limits, once is enough
                dm = [];
            end
        end

        %%-----  data model methods  -----
        function dm = data_model_update(obj, mm, nm, dm, mpopt)
            nm.pf_data_model_update(mm, nm, dm, mpopt);
        end

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

        function nm = network_model_convert_x(obj, mm, nm)
            %% convert solved state from math model to network model soln
            [nm.soln.v, nm.soln.z, nm.soln.x] = nm.pf_convert_x(mm.soln.x, ...
                                                mm.get_userdata('aux_data'));
        end

        %%-----  mathematical model methods  -----
        function mm = math_model_build_pre(obj, mm, nm, dm, mpopt)
            mm.userdata.aux_data = nm.pf_aux_data(dm, mpopt);
        end

        function obj = math_model_add_vars(obj, mm, nm, dm, mpopt)
            nm.pf_add_vars(mm, nm, dm, mpopt);
        end

        function obj = math_model_add_constraints(obj, mm, nm, dm, mpopt)
            nm.pf_add_constraints(mm, nm, dm, mpopt);
        end

        function opt = math_model_opt(obj, mm, nm, dm, mpopt)
            opt = nm.pf_solve_opts(mm, dm, mpopt);
            obj.mm_opt = opt;
        end
    end     %% methods
end         %% classdef
