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
    end     %% methods
end         %% classdef
