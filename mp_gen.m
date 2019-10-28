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

        function obj = build_params(obj, asm, mpc)
            %% define constants
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            nn = asm.getN('node');
            ng = obj.nk;
            b2n = asm.node.data.ID2idx.bus;     %% bus to node map

            %% what nodes are the gens at
            gnodes = b2n(mpc.gen(:, GEN_BUS));
            Cg = sparse(gnodes, (1:ng)', 1, nn, ng);    %% connection matrix
                                                        %% element i, j is 1 if
                                                        %% gen j is at node i
            obj.C = { Cg };
        end
    end     %% methods
end         %% classdef
