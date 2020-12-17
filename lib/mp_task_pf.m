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

%     properties
%     end
    
    methods
        %%-----  constructor  -----
        function obj = mp_task_pf()
            %% call parent constructor
            obj@mp_task();

            obj.tag = 'PF';
            obj.name = 'Power Flow';
        end

        %%-----  data model methods  -----
        function dm = data_model_update(obj, mm, nm, dm, mpopt)
            %% e.g. update data model with network model solution
            if mpopt.verbose, fprintf('-- %s data_model_update()\n', obj.tag); end
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

        function nm = network_model_update(obj, mm, nm)
            %% convert back to complex voltage vector
            ad = mm.get_userdata('power_flow_aux_data');
            [nm.soln.v, nm.soln.z] = nm.pf_convert_x(mm.soln.x, ad);
        end

        %%-----  mathematical model methods  -----
        function mm = math_model_build_pre(obj, mm, nm, dm, mpopt)
            mm.userdata.power_flow_aux_data = ...
                nm.power_flow_aux_data(dm, mpopt);
        end

        function obj = math_model_add_vars(obj, mm, nm, dm, mpopt)
            nm.add_pf_vars(nm, mm, dm, mpopt);
        end

        function obj = math_model_add_constraints(obj, mm, nm, dm, mpopt)
            nm.add_pf_constraints(nm, mm, dm, mpopt);
        end

        function opt = math_model_opt(obj, mm, nm, dm, mpopt)
            opt = nm.solve_opts_power_flow(mm, dm, mpopt);
            obj.mm_opt = opt;
        end
    end     %% methods
end         %% classdef
