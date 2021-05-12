classdef dme_gen_mpc2_node_test < dme_gen_mpc2
%DME_GEN_MPC2_NODE_TEST  MATPOWER data model gen table for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus_type    % 1 = nld, 2 = ld (all gens)
    end     %% properties

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
            dme_nld = dm.elm_by_name('bus_nld');
            dme_ld  = dm.elm_by_name('bus_ld');
            if isempty(dme_nld)
                b2i_nld = [];
            else
                b2i_nld = dme_nld.ID2i;     %% bus_nld num to idx mapping
            end
            if isempty(dme_ld)
                b2i_ld
            else
                b2i_ld = dme_ld.ID2i;       %% bus_ld num to idx mapping
            end
            n1 = length(b2i_nld);
            n2 = length(b2i_ld);
            if n1 < n2
                b2i_nld(n2, 1) = 0;
            end
            if n2 < n1
                b2i_ld(n1, 1) = 0;
            end
            b2i = b2i_nld + b2i_ld;

            %% set bus index vectors for port connectivity
            tab = obj.get_table(dm);
            k_nld = find(b2i_nld(tab(:, GEN_BUS)));
            k_ld  = find(b2i_ld( tab(:, GEN_BUS)));
            obj.bus = b2i(tab(:, GEN_BUS));
            obj.bus_type = zeros(size(obj.bus));
            obj.bus_type(k_nld) = 1;
            obj.bus_type(k_ld ) = 2;
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            dme_nld = dm.elm_by_name('bus_nld');
            dme_ld  = dm.elm_by_name('bus_ld');
            if ~isempty(dme_nld)
                bs_nld = dme_nld.status;    %% bus_nld status

                %% update status of gens at isolated/offline buses
                k_nld = find(obj.bus_type == 1);
                obj.status(k_nld) = obj.status(k_nld) & bs_nld(obj.bus(k_nld));
            end
            if ~isempty(dme_ld)
                bs_ld = dme_ld.status;      %% bus_ld status

                %% update status of gens at isolated/offline buses
                k_ld  = find(obj.bus_type == 2);
                obj.status(k_ld ) = obj.status(k_ld ) & bs_ld( obj.bus(k_ld ));
            end

            %% call parent to fill in on/off
            update_status@dme_gen(obj, dm);
        end
    end     %% methods
end         %% classdef
