classdef dmce_bus_nld_mpc2_node_test < dmce_bus_mpc2 % & dmce_bus
%DMCE_BUS_NLD_MPC2_NODE_TEST  Data model converter for bus elements for MATPOWER case v2.

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus
    end     %% properties%     end     %% properties

    methods
        function obj = dmce_bus_nld_mpc2_node_test()
            obj.name = 'bus_nld';
        end

        function [nr, nc, r] = get_size(obj, mpc)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            tab = mpc.(obj.table);
            r = find(~tab(:, PD) & ~tab(:, QD) & ~tab(:, GS) & ~tab(:, BS));
            obj.bus = r;
            nr = size(r, 1);
            nc = size(tab, 2);          %% use nc of default table
        end
    end     %% methods
end         %% classdef
