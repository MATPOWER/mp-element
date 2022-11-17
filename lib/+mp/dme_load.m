classdef dme_load < mp.dm_element
%MP.DME_LOAD  MATPOWER data model class for load data

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus     %% bus index vector (all loads)
        pd      %% active power demand (p.u.) for constant power loads that are on
        qd      %% reactive power demand (p.u.) for constant power loads that are on
        pd_i    %% active power demand (p.u.) for constant current loads that are on
        qd_i    %% reactive power demand (p.u.) for constant current loads that are on
        pd_z    %% active power demand (p.u.) for constant impedance loads that are on
        qd_z    %% reactive power demand (p.u.) for constant impedance loads that are on
    end     %% properties

    methods
        function name = name(obj)
            name = 'load';
        end

        function label = label(obj)
            label = 'Load';
        end

        function label = labels(obj)
            label = 'Loads';
        end

        function name = cxn_type(obj)
            name = 'bus';
        end

        function name = cxn_idx_prop(obj)
            name = 'bus';
        end

        function names = main_table_var_names(obj)
            names = horzcat( main_table_var_names@mp.dm_element(obj), ...
                {'bus', 'pd', 'qd', 'pd_i', 'qd_i', 'pd_z', 'qd_z', ...
                'p', 'q'});
        end

        function vars = export_vars(obj)
            vars = {};
        end

        function nr = count(obj, dm)
            nr = count@mp.dm_element(obj, dm);
            if nr
                obj.bus = obj.tab.source_uid;
            end
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elements.bus.tab.status;    %% bus status

            %% update status of loads at isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@mp.dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            obj.pd   = obj.tab.pd(obj.on) / dm.base_mva;
            obj.qd   = obj.tab.qd(obj.on) / dm.base_mva;
            obj.pd_i = obj.tab.pd_i(obj.on) / dm.base_mva;
            obj.qd_i = obj.tab.qd_i(obj.on) / dm.base_mva;
            obj.pd_z = obj.tab.pd_z(obj.on) / dm.base_mva;
            obj.qd_z = obj.tab.qd_z(obj.on) / dm.base_mva;
        end

        function TorF = pp_have_section_sum(obj, mpopt, pp_args)
            TorF = true;
        end

        function obj = pp_data_sum(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %% call parent
            pp_data_sum@mp.dm_element(obj, dm, rows, out_e, mpopt, fd, pp_args);

            %% print load summary
            fprintf(fd, '  %-29s %12.1f MW', 'Total load', ...
                                            sum(obj.tab.p));
            if mpopt.model(1) ~= 'D'    %% AC model
                fprintf(fd, ' %12.1f MVAr', sum(obj.tab.q));
            end
            fprintf(fd, '\n');
            if obj.n ~= obj.nr
                fprintf(fd, '  %-29s %12.1f MW', '  online', ...
                                                sum(obj.tab.p(obj.on)));
                if mpopt.model(1) ~= 'D'    %% AC model
                    fprintf(fd, ' %12.1f MVAr', sum(obj.tab.q(obj.on)));
                end
                fprintf(fd, '\n');
            end
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            h = [ pp_get_headers_det@mp.dm_element(obj, dm, out_e, mpopt, pp_args) ...
                {   '                             Power Consumption', ...
                    'Load ID    Bus ID   Status   P (MW)   Q (MVAr)', ...
                    '--------  --------  ------  --------  --------' } ];
            %%       1234567 123456789 -----1 12345678.0 1234567.9
        end

        function f = pp_get_footers_det(obj, dm, out_e, mpopt, pp_args)
            f = {'                            --------  --------',
                sprintf('%18s Total:%10.1f %9.1f', ...
                    '', sum(obj.tab.p(obj.on)), sum(obj.tab.q(obj.on)))};
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            if obj.tab.status(k) && abs(obj.tab.p(k)) > 1e-5
                p = sprintf('%10.1f', obj.tab.p(k));
            else
                p = '       -  ';
            end
            if obj.tab.status(k) && abs(obj.tab.q(k)) > 1e-5
                q = sprintf('%9.1f', obj.tab.q(k));
            else
                q = '      -  ';
            end
            str = sprintf('%7d %9d %6d %10s %9s', ...
                obj.tab.uid(k), obj.tab.bus(k), obj.tab.status(k), p, q);
        end
    end     %% methods
end         %% classdef
