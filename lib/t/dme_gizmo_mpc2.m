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
        function obj = update_status(obj, dm)
            %% define constants
            BUS1 = 1;
            BUS2 = 2;
            BUS3 = 3;
            
            dm_bus = dm.elm_by_name('bus');
            bs = dm_bus.status;     %% bus status
            b2i = dm_bus.ID2i;      %% bus num to idx mapping

            %% update status of gizmoes connected to isolated/offline buses
            tab = obj.get_table(dm);
            obj.status = obj.status & bs(b2i(tab(:, BUS1))) & ...
                                      bs(b2i(tab(:, BUS2))) & ...
                                      bs(b2i(tab(:, BUS3)));

            %% call parent to fill in on/off
            update_status@dme_gizmo(obj);
        end
    end     %% methods
end         %% classdef
