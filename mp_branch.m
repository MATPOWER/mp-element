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
        function obj = mp_branch(varargin)
            obj@mp_element(varargin{:});
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

        function add_opf_constraints(obj, asm, om, mpc, mpopt)
            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

            %% branch voltage angle difference limits
            nb = size(mpc.bus, 1);      %% number of buses
            [Aang, lang, uang, iang] = makeAang(mpc.baseMVA, mpc.branch, nb, mpopt);
%             nang = length(iang);
            om.add_lin_constraint('ang', Aang, lang, uang, {'Va'});
            om.userdata.iang = iang;
        end
    end     %% methods
end         %% classdef
