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
        function nm_class = network_model_class(obj, dm, mpopt)
            nm_class = obj.network_model_class_override(dm, mpopt);
            if isempty(nm_class)
                switch upper(mpopt.model)
                    case 'AC'
                        if mpopt.pf.v_cartesian
                            if mpopt.pf.current_balance
                                nm_class = @mpe_network_acci;
                            else
                                nm_class = @mpe_network_accs;
%                                nm_class = @mpe_network_accs_test_nln;
                            end
                        else
                            if mpopt.pf.current_balance
                                nm_class = @mpe_network_acpi;
                            else
                                nm_class = @mpe_network_acps;
%                                nm_class = @mpe_network_acps_test_nln;
                            end
                        end
                    case 'DC'
                        nm_class = @mpe_network_dc;
                end
            end
        end

        function nm = network_model_update(obj, mm, nm)
            %% convert back to complex voltage vector
            ad = mm.get_userdata('power_flow_aux_data');
            [nm.soln.v, nm.soln.z] = nm.pfx2vz(mm.soln.x, ad);
        end

        %%-----  mathematical model methods  -----
        function mm = math_model_create_pre(obj, mm, nm, dm, mpopt)
            ad = nm.power_flow_aux_data(dm, mpopt);
            mm.userdata.power_flow_aux_data = ad;
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
