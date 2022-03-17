classdef dme_load3p < dm_element
%DME_LOAD3P  MATPOWER data model class for load data

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus     %% bus index vector (all loads)
        pd1     %% phase 1 active power demand (p.u.) for loads that are on
        pd2     %% phase 2 active power demand (p.u.) for loads that are on
        pd3     %% phase 3 active power demand (p.u.) for loads that are on
        pf1     %% phase 1 power factor for loads that are on
        pf2     %% phase 2 power factor for loads that are on
        pf3     %% phase 3 power factor for loads that are on
    end     %% properties

    methods
        function name = name(obj)
            name = 'load3p';
        end

        function label = label(obj)
            label = '3-ph Load';
        end

        function label = labels(obj)
            label = '3-ph Loads';
        end

        function name = cxn_type(obj)
            name = 'bus3p';
        end

        function name = cxn_idx_prop(obj)
            name = 'bus';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus', 'pd1', 'pd2', 'pd3', 'pf1', 'pf2', 'pf3'});
        end

%         function vars = export_vars(obj, task)
%             vars = {};
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
            bs = dm.elements.bus3p.tab.status;  %% bus status

            %% update status of loads at isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            obj.pd1 = obj.tab.pd1(obj.on) / dm.base_kva;
            obj.pd2 = obj.tab.pd2(obj.on) / dm.base_kva;
            obj.pd3 = obj.tab.pd3(obj.on) / dm.base_kva;
            obj.pf1 = obj.tab.pf1(obj.on);
            obj.pf2 = obj.tab.pf2(obj.on);
            obj.pf3 = obj.tab.pf3(obj.on);
        end

        function obj = pp_data_sum(obj, dm, rows, out_e, mpopt, fd, varargin)
            %% call parent
            pp_data_sum@dm_element(obj, dm, rows, out_e, mpopt, fd, varargin{:});

            %% print generation summary
            fprintf(fd, '  %-29s %12.1f kW\n', 'Total 3-ph load', ...
                sum(obj.tab.pd1(obj.on)) + ...
                sum(obj.tab.pd2(obj.on)) + ...
                sum(obj.tab.pd3(obj.on)));
        end

        function TorF = pp_have_section_det(obj, mpopt, varargin)
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, varargin)
            h = {   '  3-ph      3-ph             Phase A Power     Phase B Power     Phase C Power', ...
                    'Load ID    Bus ID   Status   (kW)     (PF)     (kW)     (PF)     (kW)     (PF)', ...
                    '--------  --------  ------  -------  ------   -------  ------   -------  ------' };
            %%       1234567 123456789 -----1 1234567.90 12.4567 123456.89 12.4567 123456.89 12.4567
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, varargin)
            str = sprintf('%7d %9d %6d %10.2f %7.4f %9.2f %7.4f %9.2f %7.4f', ...
                obj.tab.uid(k), obj.tab.bus(k), obj.tab.status(k), ...
                obj.tab.pd1(k), obj.tab.pf1(k), ...
                obj.tab.pd2(k), obj.tab.pf2(k), ...
                obj.tab.pd3(k), obj.tab.pf3(k) );
        end
    end     %% methods
end         %% classdef
