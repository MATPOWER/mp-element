classdef (Abstract) nme_bus3p < mp.nm_element

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function name = name(obj)
            name = 'bus3p';
        end

        function nn = nn(obj)
            nn = 3;     %% creates 3 nodes per element
        end

        function [ref, pv, pq] = node_types(obj, nm, dm, idx)
            %% ntv = obj.node_types(nm, dm, idx)
            %% [ref, pv, pq] = obj.node_types(nm, dm, idx)
            dme = obj.data_model_element(dm);
            if nargout > 1
                ref = find(dme.type == mp.NODE_TYPE.REF);   %% ref node indices
                pv  = find(dme.type == mp.NODE_TYPE.PV);    %% PV node indices
                pq  = find(dme.type == mp.NODE_TYPE.PQ);    %% PQ node indices
            else    %% ntv - node type vector
                ref = dme.type;
            end
        end

%         function set_node_type_ref(obj, nm, dm, idx)
%             obj.data_model_element(dm).set_bus_type_ref(dm, idx);
%         end
% 
%         function set_node_type_pv(obj, nm, dm, idx)
%             obj.data_model_element(dm).set_bus_type_pv(dm, idx);
%         end
% 
%         function set_node_type_pq(obj, nm, dm, idx)
%             obj.data_model_element(dm).set_bus_type_pq(dm, idx);
%         end
    end     %% methods
end         %% classdef
