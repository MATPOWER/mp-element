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
        %% constructor
        function obj = mp_task_pf()
            %% call parent constructor
            obj@mp_task();

            obj.tag = 'PF';
            obj.name = 'Power Flow';
        end

        function nm_class = network_model_class(obj, dm, mpopt)
            switch upper(mpopt.model)
                case 'AC'
                    if mpopt.pf.v_cartesian
                        if mpopt.pf.current_balance
                            nm_class = @mpe_network_acci;
                        else
                            nm_class = @mpe_network_accs;
        %                    nm_class = @mpe_network_accs_test_nln;
                        end
                    else
                        if mpopt.pf.current_balance
                            nm_class = @mpe_network_acpi;
                        else
                            nm_class = @mpe_network_acps;
        %                    nm_class = @mpe_network_acps_test_nln;
                        end
                    end
                case 'DC'
                    nm_class = @mpe_network_dc;
            end
        end

        function ad = add_aux_data(obj, mm, nm, dm, mpopt)
            ad = nm.power_flow_aux_data(dm.mpc, mpopt);
            mm.userdata.power_flow_aux_data = ad;
        end

        function obj = add_vars(obj, mm, nm, dm, mpopt)
            nm.add_pf_vars(nm, mm, dm.mpc, mpopt);
        end

        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            nm.add_pf_constraints(nm, mm, dm.mpc, mpopt);
        end

        function opt = add_mm_opt(obj, mpopt)
            mm_opt = struct('verbose', mpopt.verbose);
            %% more to be added here

            obj.mm_opt = mm_opt;
        end

        function nm = mm2nm(obj, mm, nm)
            fprintf('-- mp_task_pf.mm2nm()\n');

            %% convert back to complex voltage vector
            x = mm.soln.x;
            out = mm.soln.output;
            ad = mm.get_userdata('power_flow_aux_data');
            [v_, z_] = nm.pfx2vz(x, ad);

            if isfield(out, 'iterations')
                i = out.iterations;
            else
                i = -1;
            end
        end

        function dm = nm2dm(obj, nm, dm, mpopt)
            fprintf('-- mp_task_pf.nm2dm()\n');
        end
    end     %% methods
end         %% classdef
