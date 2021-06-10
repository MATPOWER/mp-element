classdef dme_branch_mpc2_node_test < dme_branch_mpc2
%DME_BRANCH_MPC2_NODE_TEST  MATPOWER data model branch table for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020-2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus_elm_types = {'bus_nld', 'bus_ld'};
        nbet        %% number of bus element types
        fbus_etv    %% from bus element type vector (all gens), 1 = nld, 2 = ld
        tbus_etv    %%  to  bus element type vector (all gens), 1 = nld, 2 = ld
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
            obj.nbet = length(obj.bus_elm_types);
            for k = obj.nbet:-1:1
                if dm.elements.is_index_name(obj.bus_elm_types{k})
                    bus_dme{k} = dm.elements.(obj.bus_elm_types{k});
                    b2i_k{k} = bus_dme{k}.ID2i; %% bus element num to idx mapping
                else
                    bus_dme{k} = [];
                    b2i_k{k} = [];
                end
                n(k) = length(b2i_k{k});
            end

            %% expand individual b2i mappings to be same dimension
            n_max = max(n);
            b2i = zeros(n_max, 1);
            for k = obj.nbet:-1:1
                if n(k) < n_max
                    b2i_k{k}(n_max, 1) = 0;
                end
                b2i = b2i + b2i_k{k};
            end

            %% set bus index vectors for port connectivity
            tab = obj.get_table(dm);
            obj.fbus = b2i(tab(:, F_BUS));
            obj.tbus = b2i(tab(:, T_BUS));
            obj.fbus_etv = zeros(size(obj.fbus));
            obj.tbus_etv = zeros(size(obj.tbus));
            for k = 1:obj.nbet
                fk = find(b2i_k{k}(tab(:, F_BUS)));
                tk = find(b2i_k{k}(tab(:, T_BUS)));
                if ~isempty(bus_dme{k})
                    obj.fbus_etv(fk) = bus_dme{k}.bus_eti;
                    obj.tbus_etv(tk) = bus_dme{k}.bus_eti;
                end
            end
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            for k = 1:obj.nbet
                if dm.elements.is_index_name(obj.bus_elm_types{k})
                    bus_dme = dm.elements.(obj.bus_elm_types{k});
                    bs = bus_dme.status;    %% bus element status

                    %% update status of branches connected to isolated/offline buses
                    fk = find(obj.fbus_etv == bus_dme.bus_eti);
                    tk = find(obj.tbus_etv == bus_dme.bus_eti);
                    obj.status(fk) = obj.status(fk) & bs(obj.fbus(fk));
                    obj.status(tk) = obj.status(tk) & bs(obj.tbus(tk));
                end
            end

            %% call parent to fill in on/off
            update_status@dme_branch(obj, dm);
        end
    end     %% methods
end         %% classdef
