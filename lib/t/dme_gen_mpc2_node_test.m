classdef dme_gen_mpc2_node_test < dme_gen_mpc2
%DME_GEN_MPC2_NODE_TEST  MATPOWER data model gen table for MATPOWER case format v2

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
%         function obj = dme_gen_mpc2_node_test()
%             obj@dme_gen_mpc2();     %% call parent constructor
%         end

        function obj = initialize(obj, dm)
            initialize@dme_gen(obj, dm);    %% call parent

            %% define named indices into data matrices
            [GEN_BUS] = idx_gen;

            %% get bus mapping info
            b2i = dm.elm_by_name('bus_ld').ID2i;    %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            tab = obj.get_table(dm);
            obj.bus = b2i(tab(:, GEN_BUS));
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elm_by_name('bus_ld').status;   %% bus status

            %% update status of gens at isolated/offline buses
            obj.status = obj.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@dme_gen(obj, dm);
        end
    end     %% methods
end         %% classdef
