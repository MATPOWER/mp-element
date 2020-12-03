classdef mpe_load < mp_element

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
        function obj = mpe_load()
            obj@mp_element();
            obj.name = 'load';
            obj.dm_table = 'bus';
            obj.np = 1;             %% this is a 1 port element
        end

        function nk = count(obj, mpc)
            if isfield(mpc, obj.dm_table) && ~isempty(mpc.(obj.dm_table))
                obj.busidx = obj.load_bus(mpc);
                nk = length(obj.busidx);
                obj.nk = nk;    %% update the count stored internally
            else
                nk = 0;
            end
        end

        function obj = build_params(obj, nm, mpc)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% incidence matrices
            IDs = mpc.bus(obj.busidx, BUS_I);       %% bus IDs
            nidx = nm.node.data.ID2idx.bus(IDs);    %% node indexes
            obj.C = obj.incidence_matrix(nm.getN('node'), nidx);
            obj.D = obj.incidence_matrix(nm.getN('state'));
        end
    end     %% methods
end         %% classdef
