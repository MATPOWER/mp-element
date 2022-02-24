classdef dme_load < dm_element
%DME_LOAD  MATPOWER data model class for load data

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
        %% constructor
        function obj = dme_load()
            obj@dm_element();   %% call parent constructor
            obj.cxn_type = 'bus';
            obj.cxn_idx_prop = 'bus';
        end

        function name = name(obj)
            name = 'load';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus', 'pd', 'qd', 'pd_i', 'qd_i', 'pd_z', 'qd_z', ...
                'p', 'q'});
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
            bs = dm.elements.bus.status;    %% bus status

            %% update status of loads at isolated/offline buses
            obj.status = obj.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            obj.pd   = obj.tab.pd(obj.on) / dm.base_mva;
            obj.qd   = obj.tab.qd(obj.on) / dm.base_mva;
            obj.pd_i = obj.tab.pd_i(obj.on) / dm.base_mva;
            obj.qd_i = obj.tab.qd_i(obj.on) / dm.base_mva;
            obj.pd_z = obj.tab.pd_z(obj.on) / dm.base_mva;
            obj.qd_z = obj.tab.qd_z(obj.on) / dm.base_mva;
        end
    end     %% methods
end         %% classdef
