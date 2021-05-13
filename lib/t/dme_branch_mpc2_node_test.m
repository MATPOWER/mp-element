classdef dme_branch_mpc2_node_test < dme_branch_mpc2
%DME_BRANCH_MPC2_NODE_TEST  MATPOWER data model branch table for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        fbus_type    % 1 = nld, 2 = ld (all gens)
        tbus_type    % 1 = nld, 2 = ld (all gens)
    end     %% properties

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
            fk_nld = find(b2i_nld(tab(:, F_BUS)));
            tk_nld = find(b2i_nld(tab(:, T_BUS)));
            fk_ld  = find(b2i_ld( tab(:, F_BUS)));
            tk_ld  = find(b2i_ld( tab(:, T_BUS)));
            obj.fbus = b2i(tab(:, F_BUS));
            obj.tbus = b2i(tab(:, T_BUS));
            obj.fbus_type = zeros(size(obj.fbus));
            obj.tbus_type = zeros(size(obj.tbus));
            if ~isempty(dme_nld)
                obj.fbus_type(fk_nld) = dme_nld.bus_class;
                obj.tbus_type(tk_nld) = dme_nld.bus_class;
            end
            if ~isempty(dme_ld)
                obj.fbus_type(fk_ld ) = dme_ld.bus_class;
                obj.tbus_type(tk_ld ) = dme_ld.bus_class;
            end
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            dme_nld = dm.elm_by_name('bus_nld');
            dme_ld  = dm.elm_by_name('bus_ld');
            if ~isempty(dme_nld)
                bs_nld = dme_nld.status;    %% bus_nld status

                %% update status of branches connected to isolated/offline buses
                fk_nld = find(obj.fbus_type == dme_nld.bus_class);
                tk_nld = find(obj.tbus_type == dme_nld.bus_class);
                obj.status(fk_nld) = obj.status(fk_nld) & bs_nld(obj.fbus(fk_nld));
                obj.status(tk_nld) = obj.status(tk_nld) & bs_nld(obj.tbus(tk_nld));
            end
            if ~isempty(dme_ld)
                bs_ld = dme_ld.status;      %% bus_ld status

                %% update status of branches connected to isolated/offline buses
                fk_ld  = find(obj.fbus_type == dme_ld.bus_class);
                tk_ld  = find(obj.tbus_type == dme_ld.bus_class);
                obj.status(fk_ld ) = obj.status(fk_ld ) & bs_ld( obj.fbus(fk_ld ));
                obj.status(tk_ld ) = obj.status(tk_ld ) & bs_ld( obj.tbus(tk_ld ));
            end

            %% call parent to fill in on/off
            update_status@dme_branch(obj, dm);
        end
    end     %% methods
end         %% classdef
