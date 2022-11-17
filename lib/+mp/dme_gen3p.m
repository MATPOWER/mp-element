classdef dme_gen3p < mp.dm_element
%MP.DME_GEN3P  MATPOWER data model class for 3-phase gen data

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
        function name = name(obj)
            name = 'gen3p';
        end

        function label = label(obj)
            label = '3-ph Generator';
        end

        function label = labels(obj)
            label = '3-ph Generators';
        end

        function name = cxn_type(obj)
            name = 'bus3p';
        end

        function name = cxn_idx_prop(obj)
            name = 'bus';
        end

        function names = main_table_var_names(obj)
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'bus', 'vm1_setpoint', 'vm2_setpoint', 'vm3_setpoint', ...
                'pg1', 'pg2', 'pg3', 'qg1', 'qg2', 'qg3'});
        end

%         function vars = export_vars(obj)
%             vars = {'pg1', 'pg2', 'pg3', 'qg1', 'qg2', 'qg3'};
%         end

        function obj = initialize(obj, dm)
            initialize@mp.dm_element(obj, dm); %% call parent

            %% get bus mapping info
            b2i = dm.elements.bus3p.ID2i;   %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            obj.bus = b2i(obj.tab.bus);
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elements.bus3p.tab.status;  %% bus status

            %% update status of gens at isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            base_kva = dm.base_kva;

            gen = obj.tab;

            %% get generator parameters
            obj.pg1_start = gen.pg1(obj.on) / base_kva;
            obj.pg2_start = gen.pg2(obj.on) / base_kva;
            obj.pg3_start = gen.pg3(obj.on) / base_kva;
            obj.qg1_start = gen.qg1(obj.on) / base_kva;
            obj.qg2_start = gen.qg2(obj.on) / base_kva;
            obj.qg3_start = gen.qg3(obj.on) / base_kva;
            obj.vm1_setpoint = gen.vm1_setpoint(obj.on);
            obj.vm2_setpoint = gen.vm2_setpoint(obj.on);
            obj.vm3_setpoint = gen.vm3_setpoint(obj.on);
        end

        function TorF = pp_have_section_sum(obj, mpopt, pp_args)
            TorF = true;
        end

        function obj = pp_data_sum(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %% call parent
            pp_data_sum@mp.dm_element(obj, dm, rows, out_e, mpopt, fd, pp_args);

            %% print generation summary
            fprintf(fd, '  %-29s %12.1f kW %12.1f kVAr\n', 'Total 3-ph generation', ...
                sum(obj.tab.pg1(obj.on)) + ...
                sum(obj.tab.pg2(obj.on)) + ...
                sum(obj.tab.pg3(obj.on)), ...
                sum(obj.tab.qg1(obj.on)) + ...
                sum(obj.tab.qg2(obj.on)) + ...
                sum(obj.tab.qg3(obj.on)) );
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            h = [ pp_get_headers_det@mp.dm_element(obj, dm, out_e, mpopt, pp_args) ...
                {   '  3-ph      3-ph             Phase A Power     Phase B Power     Phase C Power', ...
                    ' Gen ID    Bus ID   Status   (kW)    (KVAr)    (kW)    (kVAr)    (kW)    (kVAr)', ...
                    '--------  --------  ------  -------  ------   -------  ------   -------  ------' } ];
            %%       1234567 123456789 -----1 1234567.90 1234.67 123456.89 1234.67 123456.89 1234.67
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            str = sprintf('%7d %9d %6d %10.2f %7.2f %9.2f %7.2f %9.2f %7.2f', ...
                obj.tab.uid(k), obj.tab.bus(k), obj.tab.status(k), ...
                obj.tab.pg1(k), obj.tab.qg1(k), ...
                obj.tab.pg2(k), obj.tab.qg2(k), ...
                obj.tab.pg3(k), obj.tab.qg3(k) );
        end
    end     %% methods
end         %% classdef
