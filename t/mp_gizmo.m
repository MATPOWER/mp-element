classdef mp_gizmo < mp_element

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
        function obj = mp_gizmo(varargin)
            obj@mp_element(varargin{:});
            obj.name = 'gizmo';
            obj.mpc_field = 'gizmo';
            obj.np = 3;             %% this is a 3 port element
            obj.nz = 2;             %% each with 2 state variables
        end

        function obj = add_states(obj, asm, mpc)
            if obj.nz > 1
                asm.init_indexed_name('state', obj.name, {obj.nz});
                for k = 1:obj.nz
                    asm.add_state(obj.name, {k}, obj.nk);
                end
            elseif obj.nz == 1
                asm.add_state(obj.name, obj.nk);
            end
        end

        function obj = build_params(obj, asm, mpc)
            gizmo = mpc.gizmo;
            nk = obj.nk;

            %% incidence matrices
            nn = asm.getN('node');
            ID2idx = asm.node.data.ID2idx.bus;
            idx1 = ID2idx(mpc.gizmo(:, 1)); %% port 1 node indexes
            idx2 = ID2idx(mpc.gizmo(:, 2)); %% port 2 node indexes
            idx3 = ID2idx(mpc.gizmo(:, 3)); %% port 3 node indexes
            obj.C = { sparse(idx1, 1:nk, 1, nn, nk), ...
                      sparse(idx2, 1:nk, 1, nn, nk), ...
                      sparse(idx3, 1:nk, 1, nn, nk) };

            nz = asm.getN('state');
            ss = asm.get_idx('state');
            obj.D = { sparse(ss.i1.(obj.name)(1):ss.iN.(obj.name)(1), 1:nk, 1, nz, nk), ...
                      sparse(ss.i1.(obj.name)(2):ss.iN.(obj.name)(2), 1:nk, 1, nz, nk) };
        end
    end     %% methods
end         %% classdef
