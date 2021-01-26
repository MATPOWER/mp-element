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

        function ntv = node_types(obj, nm, dm)
            dme = obj.data_model_element(dm);
            ntv = dme.bus_types(dm);
            if nargout > 1
                ref = dm.node_type_ref(ntv);    %% reference node indices
                pv  = dm.node_type_pv(ntv);     %% PV node indices
                pq  = dm.node_type_pq(ntv);     %% PQ node indices
            else
                ref = ntv;
            end
        end
        function set_node_type_ref(obj, nm, dm, idx)
            dme = obj.data_model_element(dm);
            dme.isref(idx) = 1;
            dme.ispv(idx)  = 0;
            dme.ispq(idx)  = 0;
        end

        function set_node_type_pq(obj, nm, dm, idx)
            dme = obj.data_model_element(dm);
            dme.isref(idx) = 0;
            dme.ispv(idx)  = 0;
            dme.ispq(idx)  = 1;
        end

    end     %% methods
end         %% classdef
