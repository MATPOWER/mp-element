classdef nme_gizmo < nm_element

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'gizmo';
%     end
    
    methods
        %% constructor
        function obj = nme_gizmo()
            obj@nm_element();
            obj.name = 'gizmo';
            obj.np = 3;             %% this is a 3 port element
            obj.nz = 2;             %% each with 2 state variables
        end

        function obj = add_states(obj, nm, dm)
            if obj.nz > 1
                nm.init_indexed_name('state', obj.name, {obj.nz});
                for k = 1:obj.nz
                    nm.add_state(obj.name, {k}, obj.nk);
                end
            elseif obj.nz == 1
                nm.add_state(obj.name, obj.nk);
            end
        end

        function obj = build_params(obj, nm, dm)
            dme = obj.data_model_element(dm);

            %% incidence matrices
            nidx = nm.get_node_idx('bus');      %% node indices for 'bus'
            sidx = nm.get_state_idx(obj.name);  %% state indices for 'gizmo'
            idx1 = nidx(dme.bus1);              %% gizmo port 1 node indices
            idx2 = nidx(dme.bus2);              %% gizmo port 2 node indices
            idx3 = nidx(dme.bus3);              %% gizmo port 3 node indices
            obj.C = obj.incidence_matrix(nm.getN('node'), idx1, idx2, idx3);
            obj.D = obj.incidence_matrix(nm.getN('state'), sidx{1}, sidx{2});
        end
    end     %% methods
end         %% classdef
