classdef dme_branch_mpc2_node_test < dme_branch_mpc2
%DME_BRANCH_MPC2_NODE_TEST  MATPOWER data model branch table for MATPOWER case format v2

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
%         %% constructor
%         function obj = dme_branch_mpc2_node_test()
%             obj@dme_branch();   %% call parent constructor
%         end

        function obj = initialize(obj, dm)
            initialize@dme_branch(obj, dm);     %% call parent

            %% define named indices into data matrices
            [F_BUS, T_BUS] = idx_brch;

            %% get bus mapping info
            b2i = dm.elm_by_name('bus_ld').ID2i;   %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            tab = obj.get_table(dm);
            obj.fbus = b2i(tab(:, F_BUS));
            obj.tbus = b2i(tab(:, T_BUS));
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elm_by_name('bus_ld').status;   %% bus status

            %% update status of branches connected to isolated/offline buses
            obj.status = obj.status & bs(obj.fbus) & ...
                                      bs(obj.tbus);

            %% call parent to fill in on/off
            update_status@dme_branch(obj, dm);
        end
    end     %% methods
end         %% classdef
