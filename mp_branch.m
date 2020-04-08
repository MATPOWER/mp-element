classdef mp_branch < mp_element

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'branch';
%     end
    
    methods
        %% constructor
        function obj = mp_branch()
            obj@mp_element();
            obj.name = 'branch';
            obj.mpc_field = 'branch';
            obj.np = 2;             %% this is a 2 port element
        end

        function obj = build_params(obj, asm, mpc)
            %% define constants
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS] = idx_brch;

            branch = mpc.branch;
            nl = obj.nk;

            %% incidence matrices
            nn = asm.getN('node');
            nl = obj.nk;

            fIDs = mpc.branch(:, F_BUS);            %% "from" bus IDs
            tIDs = mpc.branch(:, T_BUS);            %% "to" bus IDs
            fidx = asm.node.data.ID2idx.bus(fIDs);  %% "from" node indexes
            tidx = asm.node.data.ID2idx.bus(tIDs);  %% "to" node indexes

            obj.C = { sparse(fidx, 1:nl, 1, nn, nl), ...
                      sparse(tidx, 1:nl, 1, nn, nl) };
        end
    end     %% methods
end         %% classdef
