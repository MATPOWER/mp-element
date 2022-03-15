classdef dme_branch < dm_element
%DME_BRANCH  MATPOWER data model class for branch data

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        fbus    %% bus index vector for "from" port (port 1) (all branches)
        tbus    %% bus index vector for "to" port (port 2) (all branches)
        r       %% series resistance (p.u.) for branches that are on
        x       %% series reactance (p.u.) for branches that are on
        g_fr    %% shunt conductance (p.u.) at "from" end for branches that are on
        g_to    %% shunt conductance (p.u.) at "to" end for branches that are on
        b_fr    %% shunt susceptance (p.u.) at "from" end for branches that are on
        b_to    %% shunt susceptance (p.u.) at "to" end for branches that are on
        tm      %% transformer off-nominal turns ratio for branches that are on
        ta      %% xformer phase-shift angle (radians) for branches that are on
        rate_a  %% long term flow limit (p.u.) for branches that are on
    end     %% properties

    methods
        function name = name(obj)
            name = 'branch';
        end

        function label = label(obj)
            label = 'Branch';
        end

        function label = labels(obj)
            label = 'Branches';
        end

        function name = cxn_type(obj)
            name = 'bus';
        end

        function name = cxn_idx_prop(obj)
            name = {'fbus', 'tbus'};
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus_fr', 'bus_to', 'r', 'x', 'g_fr', 'b_fr', ...
                'g_to', 'b_to', 'sm_ub_a', 'sm_ub_b', 'sm_ub_c', ...
                'cm_ub_a', 'cm_ub_b', 'cm_ub_c', 'vad_lb', 'vad_ub', ...
                'tm', 'ta', ... %% remove these when we separate out xformers
                'pl_fr', 'ql_fr', 'pl_to', 'ql_to'} );
        end

        function vars = export_vars(obj, task)
            vars = {'ql_to', 'pl_to', 'ql_fr', 'pl_fr'};
        end

        function obj = initialize(obj, dm)
            initialize@dm_element(obj, dm);     %% call parent

            %% get bus mapping info
            b2i = dm.elements.bus.ID2i;     %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            obj.fbus = b2i(obj.tab.bus_fr);
            obj.tbus = b2i(obj.tab.bus_to);
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elements.bus.tab.status;    %% bus status

            %% update status of branches connected to isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.fbus) & ...
                                              bs(obj.tbus);

            %% call parent to fill in on/off
            update_status@dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            obj.r  = obj.tab.r(obj.on);
            obj.x  = obj.tab.x(obj.on);
            obj.g_fr  = obj.tab.g_fr(obj.on);
            obj.b_fr  = obj.tab.b_fr(obj.on);
            obj.g_to  = obj.tab.g_to(obj.on);
            obj.b_to  = obj.tab.b_to(obj.on);
            obj.tm = obj.tab.tm(obj.on);
            obj.ta = obj.tab.ta(obj.on) * pi/180;
            obj.rate_a = obj.tab.sm_ub_a(obj.on) / dm.base_mva;
        end
    end     %% methods
end         %% classdef
