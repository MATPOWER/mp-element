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

        function obj = build_params(obj, nm, dm)
            dme = obj.data_model_element(dm);

            %% incidence matrices
            nidxs = obj.node_indices(nm, dm, 'bus', {'bus1', 'bus2', 'bus3'});
            sidx = nm.get_state_idx(obj.name);  %% state indices for 'gizmo'
            obj.C = obj.incidence_matrix(nm.getN('node'), nidxs{:});
            obj.D = obj.incidence_matrix(nm.getN('state'), sidx{1}, sidx{2});
        end
    end     %% methods
end         %% classdef
