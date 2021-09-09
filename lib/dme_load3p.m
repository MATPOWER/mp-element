classdef dme_load3p < dm_element
%DME_LOAD3P  MATPOWER data model class for load data

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
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
        %% constructor
        function obj = dme_load3p()
            obj@dm_element();   %% call parent constructor
            obj.name = 'load3p';
            obj.cxn_type = 'bus3p';
            obj.cxn_idx_prop = 'bus';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus', 'pd1', 'pd2', 'pd3', 'pf1', 'pf2', 'pf3'});
        end

%         function vars = export_vars(obj, task)
%             switch task
%                 case 'PF'
%                     vars = {};
%                 case 'CPF'
%                     vars = {'pd', 'qd', 'pd_i', 'qd_i', 'pd_z', 'qd_z'};
%                 case 'OPF'
%                     vars = {};
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
            obj.pd1 = obj.tab.pd1(obj.on) / dm.base_kva;
            obj.pd2 = obj.tab.pd2(obj.on) / dm.base_kva;
            obj.pd3 = obj.tab.pd3(obj.on) / dm.base_kva;
            obj.pf1 = obj.tab.pf1(obj.on);
            obj.pf2 = obj.tab.pf2(obj.on);
            obj.pf3 = obj.tab.pf3(obj.on);
        end
    end     %% methods
end         %% classdef
