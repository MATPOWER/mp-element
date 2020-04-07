classdef mp_shunt < mp_element

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        busidx = [];
    end

    methods
        %% constructor
        function obj = mp_shunt()
            obj@mp_element();
            obj.name = 'shunt';
            obj.mpc_field = 'bus';
            obj.np = 1;             %% this is a 1 port element
        end

        function nk = count(obj, mpc)
            if isfield(mpc, obj.mpc_field) && ~isempty(mpc.(obj.mpc_field))
                obj.busidx = obj.shunt_bus(mpc);
                nk = length(obj.busidx);
                obj.nk = nk;    %% update the count stored internally
            else
                nk = 0;
            end
        end

        function obj = build_params(obj, asm, mpc)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% incidence matrices
            nn = asm.getN('node');
            nsh = obj.nk;
            k = obj.busidx;
            IDs = mpc.bus(k, BUS_I);                %% bus IDs
            nidx = asm.node.data.ID2idx.bus(IDs);   %% node indexes
            obj.C = { sparse(nidx, 1:nsh, 1, nn, nsh) };
        end
    end     %% methods
end         %% classdef
