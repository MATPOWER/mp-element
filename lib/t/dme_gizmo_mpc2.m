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

            %% get bus mapping info
            b2i = dm.elements.bus.ID2i;     %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            tab = obj.get_table(dm);
            obj.bus1 = b2i(tab(:, BUS1));
            obj.bus2 = b2i(tab(:, BUS2));
            obj.bus3 = b2i(tab(:, BUS3));
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elements.bus.status;    %% bus status

            %% update status of gizmoes connected to isolated/offline buses
            obj.status = obj.status & bs(obj.bus1) & ...
                                      bs(obj.bus2) & ...
                                      bs(obj.bus3);

            %% call parent to fill in on/off
            update_status@dme_gizmo(obj, dm);
        end
    end     %% methods
end         %% classdef
