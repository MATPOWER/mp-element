classdef dme_xfmr3p < dm_element
%DME_XFMR3P  MATPOWER data model class for 3-phase transformer data

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
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
        base_kva
        base_kv
%         g_to    %% shunt conductance (p.u.) at "to" end for branches that are on
%         b_fr    %% shunt susceptance (p.u.) at "from" end for branches that are on
%         b_to    %% shunt susceptance (p.u.) at "to" end for branches that are on
%         tm      %% transformer off-nominal turns ratio for branches that are on
%         ta      %% xformer phase-shift angle (radians) for branches that are on
%         rate_a  %% long term flow limit (p.u.) for branches that are on
    end     %% properties

    methods
        function name = name(obj)
            name = 'xfmr3p';
        end

        function name = cxn_type(obj)
            name = 'bus3p';
        end

        function name = cxn_idx_prop(obj)
            name = {'fbus', 'tbus'};
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus_fr', 'bus_to', 'r', 'x', 'base_kva', 'base_kv', ...
                 'pl1_fr', 'pf1_fr', 'pl2_fr', 'pf2_fr', 'pl3_fr', 'pf3_fr', ...
                 'pl1_to', 'pf1_to', 'pl2_to', 'pf2_to', 'pl3_to', 'pf3_to' ...
                 });
        end

%         function vars = export_vars(obj, task)
%             switch task
%                 case 'PF'
%                     vars = {'ql_to', 'pl_to', 'ql_fr', 'pl_fr'};
%                 case 'CPF'
%                     vars = {'ql_to', 'pl_to', 'ql_fr', 'pl_fr'};
%                 case 'OPF'
%                     vars = {'ql_to', 'pl_to', 'ql_fr', 'pl_fr', ...
%                         'mu_flow_fr_ub', 'mu_flow_to_ub', ...
%                         'mu_vad_lb', 'mu_vad_ub'};
%                 otherwise
%                     vars = 'all';
%             end
%         end

        function obj = initialize(obj, dm)
            initialize@dm_element(obj, dm); %% call parent

            %% get bus mapping info
            b2i = dm.elements.bus3p.ID2i;   %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            obj.fbus = b2i(obj.tab.bus_fr);
            obj.tbus = b2i(obj.tab.bus_to);
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elements.bus3p.status;  %% bus status

            %% update status of branches connected to isolated/offline buses
            obj.status = obj.status & bs(obj.fbus) & ...
                                      bs(obj.tbus);

            %% call parent to fill in on/off
            update_status@dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            obj.r = obj.tab.r(obj.on);
            obj.x = obj.tab.x(obj.on);
            obj.base_kva = obj.tab.base_kva(obj.on);
            obj.base_kv  = obj.tab.base_kv(obj.on);
        end
    end     %% methods
end         %% classdef
