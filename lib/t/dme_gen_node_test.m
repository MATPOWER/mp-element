classdef dme_gen_node_test < dme_gen
%DME_GEN_NODE_TEST  MATPOWER data model gen table for MATPOWER case format v2

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
        bus_etv    %% bus element type vector (all gens), 1 = nld, 2 = ld
    end     %% properties

    methods
        %% constructor
%         function obj = dme_gen_node_test()
%             obj@dme_gen();     %% call parent constructor
%         end

        function obj = initialize(obj, dm)
            initialize@dm_element(obj, dm);     %% call parent

            %% get bus mapping info
            obj.nbet = length(obj.bus_elm_types);
            for k = obj.nbet:-1:1
                if dm.elements.is_index_name(obj.bus_elm_types{k})
                    bus_dme{k} = dm.elements.(obj.bus_elm_types{k});
                    b2i_k{k} = bus_dme{k}.ID2i;
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
            obj.bus = b2i(obj.tab.bus);
            obj.bus_etv = zeros(size(obj.bus));
            for k = 1:obj.nbet
                gk = find(b2i_k{k}(obj.tab.bus));
                if ~isempty(bus_dme{k})
                    obj.bus_etv(gk) = bus_dme{k}.bus_eti;
                end
            end
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            for k = 1:obj.nbet
                if dm.elements.is_index_name(obj.bus_elm_types{k})
                    bus_dme = dm.elements.(obj.bus_elm_types{k});
                    bs = bus_dme.status;    %% bus element status

                    %% update status of gens at isolated/offline buses
                    gk = find(obj.bus_etv == bus_dme.bus_eti);
                    obj.status(gk) = obj.status(gk) & bs(obj.bus(gk));
                end
            end

            %% call parent to fill in on/off
            update_status@dm_element(obj, dm);
        end
    end     %% methods
end         %% classdef
