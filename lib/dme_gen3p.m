classdef dme_gen3p < dm_element
%DME_GEN3P  MATPOWER data model class for 3-phase gen data

%   MATPOWER
%   Copyright (c) 1996-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus         %% bus index vector (all gens)
        pg1_start   %% initial phase 1 active power (p.u.) for gens that are on
        pg2_start   %% initial phase 2 active power (p.u.) for gens that are on
        pg3_start   %% initial phase 3 active power (p.u.) for gens that are on
        qg1_start   %% initial phase 1 reactive power (p.u.) for gens that are on
        qg2_start   %% initial phase 2 reactive power (p.u.) for gens that are on
        qg3_start   %% initial phase 3 reactive power (p.u.) for gens that are on
        vm1_setpoint%% phase 1 generator voltage setpoint for gens that are on
        vm2_setpoint%% phase 2 generator voltage setpoint for gens that are on
        vm3_setpoint%% phase 3 generator voltage setpoint for gens that are on
    end     %% properties

    methods
        %% constructor
        function obj = dme_gen3p()
            obj@dm_element();   %% call parent constructor
            obj.cxn_type = 'bus3p';
            obj.cxn_idx_prop = 'bus';
        end

        function name = name(obj)
            name = 'gen3p';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus', 'vm1_setpoint', 'vm2_setpoint', 'vm3_setpoint', ...
                'pg1', 'pg2', 'pg3', 'pf1', 'pf2', 'pf3'});
        end

%         function vars = export_vars(obj, task)
%             switch task
%                 case 'PF'
%                     vars = {'pg', 'qg'};
%                 case 'CPF'
%                     vars = {'pg', 'qg'};
%                 case 'OPF'
%                     vars = {'vm_setpoint', 'pg', 'qg', 'mu_pg_lb', 'mu_pg_ub', 'mu_qg_lb', 'mu_qg_ub'};
%                 otherwise
%                     vars = 'all';
%             end
%         end

        function obj = initialize(obj, dm)
            initialize@dm_element(obj, dm); %% call parent

            %% get bus mapping info
            b2i = dm.elements.bus3p.ID2i;   %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            obj.bus = b2i(obj.tab.bus);
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elements.bus3p.status;  %% bus status

            %% update status of gens at isolated/offline buses
            obj.status = obj.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            base_kva = dm.base_kva;

            gen = obj.tab;

            %% get generator parameters
            obj.pg1_start = gen.pg1(obj.on) / base_kva;
            obj.pg2_start = gen.pg2(obj.on) / base_kva;
            obj.pg3_start = gen.pg3(obj.on) / base_kva;
            obj.qg1_start = obj.pg1_start .* tan(acos(gen.pf1(obj.on)));
            obj.qg2_start = obj.pg2_start .* tan(acos(gen.pf2(obj.on)));
            obj.qg3_start = obj.pg3_start .* tan(acos(gen.pf3(obj.on)));
            obj.vm1_setpoint = gen.vm1_setpoint(obj.on);
            obj.vm2_setpoint = gen.vm2_setpoint(obj.on);
            obj.vm3_setpoint = gen.vm3_setpoint(obj.on);
        end
    end     %% methods
end         %% classdef
