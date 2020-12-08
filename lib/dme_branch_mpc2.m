classdef dme_branch_mpc2 < dme_branch & dm_format_mpc2
%DME_BRANCH_MPC2  MATPOWER data model branch table for MATPOWER case format v2

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
        %% constructor
        function obj = dme_branch_mpc2()
            obj@dme_branch();   %% call parent constructor

            %% define constants
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                    TAP, SHIFT, BR_STATUS] = idx_brch;
            obj.st_col = BR_STATUS;
        end

        function obj = initialize(obj, dm)
            initialize@dme_branch(obj, dm);     %% call parent

            %% define constants
            [F_BUS, T_BUS] = idx_brch;

            %% set busIDs for connectivity
            tab = obj.get_table(dm);
            obj.fbusID = tab(:, F_BUS);
            obj.tbusID = tab(:, T_BUS);
        end

        function obj = update_status(obj, dm)
            %% get bus status/mapping info
            dme_bus = dm.elm_by_name('bus');
            bs = dme_bus.status;    %% bus status
            b2i = dme_bus.ID2i;     %% bus num to idx mapping

            %% update status of branches connected to isolated/offline buses
            obj.status = obj.status & bs(b2i(obj.fbusID)) & ...
                                      bs(b2i(obj.tbusID));

            %% call parent to fill in on/off
            update_status@dme_branch(obj, dm);
        end
    end     %% methods
end         %% classdef
