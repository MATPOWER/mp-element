classdef dme_gizmo_mpc2 < dme_gizmo & dm_format_mpc2
%DME_GIZMO_MPC2  MATPOWER data model gizmo table for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = initialize(obj, dm)
            initialize@dme_gizmo(obj, dm);      %% call parent

            %% define constants
            BUS1 = 1;
            BUS2 = 2;
            BUS3 = 3;

            %% set busIDs for connectivity
            tab = obj.get_table(dm);
            obj.bus1ID = tab(:, BUS1);
            obj.bus2ID = tab(:, BUS2);
            obj.bus3ID = tab(:, BUS3);
        end

        function obj = update_status(obj, dm)
            %% get bus status/mapping info
            dme_bus = dm.elm_by_name('bus');
            bs = dme_bus.status;    %% bus status
            b2i = dme_bus.ID2i;     %% bus num to idx mapping

            %% update status of gizmoes connected to isolated/offline buses
            obj.status = obj.status & bs(b2i(obj.bus1ID)) & ...
                                      bs(b2i(obj.bus2ID)) & ...
                                      bs(b2i(obj.bus3ID));

            %% call parent to fill in on/off
            update_status@dme_gizmo(obj, dm);
        end
    end     %% methods
end         %% classdef
