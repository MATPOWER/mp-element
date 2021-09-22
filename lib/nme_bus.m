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
            obj.nn = 1;         %% creates 1 node per element
        end

        function [ref, pv, pq] = node_types(obj, nm, dm, idx)
            %% ntv = obj.node_types(nm, dm, idx)
            %% [ref, pv, pq] = obj.node_types(nm, dm, idx)
            dme = obj.data_model_element(dm);
            if nargout > 1
                ref = find(dme.type == NODE_TYPE.REF);  %% ref node indices
                pv  = find(dme.type == NODE_TYPE.PV);   %% PV node indices
                pq  = find(dme.type == NODE_TYPE.PQ);   %% PQ node indices
            else    %% ntv - node type vector
                ref = dme.type;
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

        %%-----  OPF methods  -----
        function vm = opf_interior_vm(obj, mm, nm, dm)
            %% return vm equal to avg of clipped limits
            dme = obj.data_model_element(dm);
            vm_ub = min(dme.vm_ub, 1.5);
            vm_lb = max(dme.vm_lb, 0.5);
            vm = (vm_ub + vm_lb) / 2;
        end
    end     %% methods
end         %% classdef
