classdef dme_shunt < dm_element
%DME_SHUNT  MATPOWER data model class for shunt data

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus     %% bus index vector (all shunts)
        gs      %% shunt conductance (p.u. active power demanded at
                %% V = 1.0 p.u.) for shunts that are on
        bs      %% shunt susceptance (p.u. reactive power injected at
                %% V = 1.0 p.u.) for shunts that are on
    end     %% properties

    methods
        function name = name(obj)
            name = 'shunt';
        end

        function label = label(obj)
            label = 'Fixed Shunt';
        end

        function label = labels(obj)
            label = 'Fixed Shunts';
        end

        function name = cxn_type(obj)
            name = 'bus';
        end

        function name = cxn_idx_prop(obj)
            name = 'bus';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus', 'gs', 'bs', 'p', 'q'});
        end

        function vars = export_vars(obj, task)
            vars = {};
        end

        function nr = count(obj, dm)
            nr = count@dm_element(obj, dm);
            if nr
                obj.bus = obj.tab.source_uid;
            end
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elements.bus.tab.status;    %% bus status

            %% update status of shunts at isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            obj.gs = obj.tab.gs(obj.on) / dm.base_mva;
            obj.bs = obj.tab.bs(obj.on) / dm.base_mva;
        end

        function TorF = pp_have_section_sum(obj, mpopt, pp_args)
            TorF = true;
        end

        function obj = pp_data_sum(obj, dm, rows, out_e, mpopt, fd, pp_args)
            %% call parent
            pp_data_sum@dm_element(obj, dm, rows, out_e, mpopt, fd, pp_args);

            %% print shunt summary
            fprintf(fd, '  %-29s %12.1f MW %12.1f MVAr\n', 'Total shunt', ...
                sum(obj.tab.p), sum(obj.tab.q));
            if obj.n ~= obj.nr
                fprintf(fd, '  %-29s %12.1f MW %12.1f MVAr\n', '  online', ...
                    sum(obj.tab.p(obj.on)), sum(obj.tab.q(obj.on)));
            end
        end

        function TorF = pp_have_section_det(obj, mpopt, pp_args)
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, pp_args)
            h = [ pp_get_headers_det@dm_element(obj, dm, out_e, mpopt, pp_args) ...
                {   '                             Power Consumption', ...
                    'Shunt ID   Bus ID   Status   P (MW)   Q (MVAr)', ...
                    '--------  --------  ------  --------  --------' } ];
            %%       1234567 123456789 -----1 1234567.90 123456.89
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, pp_args)
            str = sprintf('%7d %9d %6d %10.2f %9.2f', ...
                obj.tab.uid(k), obj.tab.bus(k), obj.tab.status(k), ...
                obj.tab.p(k), obj.tab.q(k));
        end
    end     %% methods
end         %% classdef
