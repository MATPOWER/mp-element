classdef dmce_bus_ld_mpc2_node_test < dmce_bus_nld_mpc2_node_test % & dmce_bus
%DMCE_BUS_LD_MPC2_NODE_TEST  Data model converter for bus elements for MATPOWER case v2.

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = dmce_bus_ld_mpc2_node_test()
            obj.name = 'bus_ld';
        end

        function [nr, nc, r] = get_size(obj, mpc)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            tab = mpc.(obj.table);
            r = find(tab(:, PD) | tab(:, QD) | tab(:, GS) | tab(:, BS));
            obj.bus = r;
            nr = size(r, 1);
            nc = size(tab, 2);          %% use nc of default table
        end

        function vmap = table_var_map(obj, var_names, mpc)
            vmap = table_var_map@dmce_bus_nld_mpc2_node_test(obj, var_names);

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% map arguments for each name
            vmap.pd.args = PD;
            vmap.qd.args = QD;
            vmap.gs.args = GS;
            vmap.bs.args = BS;

            %% map type for each name (default mapping is -1)
            vmap.source_uid.type = 4;    %% index in mpc.bus
        end
    end     %% methods
end         %% classdef
