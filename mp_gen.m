classdef mp_gen < mp_element

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'gen';
%     end
    
    methods
        %% constructor
        function obj = mp_gen(varargin)
            obj@mp_element(varargin{:});
            obj.name = 'gen';
            obj.mpc_field = 'gen';
            obj.np = 1;             %% this is a 1 port element
            obj.nz = 1;
        end

        function obj = add_states(obj, asm, mpc)
            %% define constants
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS] = idx_gen;
    
            ng = obj.nk;            %% number of gens
            asm.add_state(obj.name, ng);
        end

        function obj = build_params(obj, asm, mpc)
            %% define constants
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS] = idx_gen;

            %% incidence matrices
            nn = asm.getN('node');
            ng = obj.nk;
            IDs = mpc.gen(:, GEN_BUS);              %% bus IDs
            nidx = asm.node.data.ID2idx.bus(IDs);   %% node indexes
            obj.C = { sparse(nidx, 1:ng, 1, nn, ng) };

            nz = asm.getN('state');
            sidx = asm.state.data.ID2idx.(obj.name);    %% state indexes
            obj.D = { sparse(sidx, 1:ng, 1, nz, ng) };
        end
    end     %% methods
end         %% classdef
