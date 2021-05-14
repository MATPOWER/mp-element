classdef nme_bus < nm_element

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'bus';
%     end

    methods
        %% constructor
        function obj = nme_bus()
            obj@nm_element();
            obj.name = 'bus';
            obj.np = 0;             %% this is a 0 port element
        end

        function obj = add_nodes(obj, nm, dm)
            dme = obj.data_model_element(dm);
            nm.add_node(obj.name, obj.nk);
        end

        function [ref, pv, pq] = node_types(obj, nm, dm)
            %% ntv = obj.node_types(nm, dm)
            %% [ref, pv, pq] = obj.node_types(nm, dm)
            dme = obj.data_model_element(dm);
            if nargout > 1
                ref = find(dme.isref);      %% reference node indices
                pv  = find(dme.ispv);       %% PV node indices
                pq  = find(dme.ispq);       %% PQ node indices
            else    %% ntv
                ref = dm.node_type_vector(dme.isref, dme.ispv, dme.ispq);
            end
        end

        function set_node_type_ref(obj, nm, dm, idx)
            obj.data_model_element(dm).set_bus_type_ref(dm, idx);
        end

        function set_node_type_pv(obj, nm, dm, idx)
            obj.data_model_element(dm).set_bus_type_pv(dm, idx);
        end

        function set_node_type_pq(obj, nm, dm, idx)
            obj.data_model_element(dm).set_bus_type_pq(dm, idx);
        end
    end     %% methods
end         %% classdef
