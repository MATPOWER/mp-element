classdef mpe_branch < mp_element

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
        function obj = mpe_branch()
            obj@mp_element();
            obj.name = 'branch';
            obj.mpc_field = 'branch';
            obj.np = 2;             %% this is a 2 port element
        end

        function obj = build_params(obj, nm, mpc)
            %% define constants
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS] = idx_brch;

            %% incidence matrices
            fIDs = mpc.branch(:, F_BUS);            %% "from" bus IDs
            tIDs = mpc.branch(:, T_BUS);            %% "to" bus IDs
            fidx = nm.node.data.ID2idx.bus(fIDs);   %% "from" node indexes
            tidx = nm.node.data.ID2idx.bus(tIDs);   %% "to" node indexes
            obj.C = obj.incidence_matrix(nm.getN('node'), fidx, tidx);
            obj.D = obj.incidence_matrix(nm.getN('state'));
        end
    end     %% methods
end         %% classdef
