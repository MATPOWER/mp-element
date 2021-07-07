classdef dmce_shunt_mpc2 < dmc_element_mpc2 % & dmce_shunt
%DMCE_SHUNT_MPC2  Data model converter for shunt elements for MATPOWER case v2.

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus
    end     %% properties

    methods
        function obj = dmce_shunt_mpc2()
            obj.name = 'shunt';
            obj.table = 'bus';
        end

        function [nr, nc, r] = get_import_size(obj, mpc)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            tab = mpc.(obj.table);
            r = find(tab(:, GS) | tab(:, BS));
            obj.bus = r;
            nr = size(r, 1);
            nc = size(tab, 2);          %% use nc of default table
        end

        function [nr, nc, r] = get_export_size(obj, dme)
            [nr, nc] = size(dme.tab);   %% use size of default table
            r = dme.tab.source_uid;     %% rows in bus matrix
        end

        function vmap = table_var_map(obj, var_names, mpc)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names, mpc);

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% map type for each name (default mapping is -1)
            vmap.uid.type        = 3;    %% consecutive IDs, starting at 1
            vmap.name.type       = 2;    %% empty char
            vmap.source_uid.type = 4;    %% index in mpc.bus
            vmap.status.type     = 1;    %% ones
            vmap.p.type          = 0;    %% zeros
            vmap.q.type          = 0;    %% zeros

            %% map arguments for each name
           %vmap.uid.args           = [];
           %vmap.name.args          = [];
           %vmap.status.args        = [];
           %vmap.source_uid.args    = [];
            vmap.bus.args           = BUS_I;
            vmap.gs.args            = GS;
            vmap.bs.args            = BS;
           %vmap.p.args             = [];
           %vmap.q.args             = [];
        end
    end     %% methods
end         %% classdef
