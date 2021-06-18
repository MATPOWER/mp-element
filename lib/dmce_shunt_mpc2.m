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
        end

        function vmap = table_var_map(obj, var_names)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names);

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% map column indices for each name (0 for exceptions)
           %vmap.idx.uid        = 0;
           %vmap.idx.name       = 0;
           %vmap.idx.status     = 0;
           %vmap.idx.source_uid = 0;
            vmap.idx.bus        = BUS_I;
            vmap.idx.gs         = GS;
            vmap.idx.bs         = BS;
           %vmap.idx.p          = 0;
           %vmap.idx.q          = 0;

            %% map table for each name (0 for default mapping)
            vmap.tab.uid        = 2;    %% consecutive IDs, starting at 1
            vmap.tab.name       = 1;    %% empty char
            vmap.tab.source_uid = 4;    %% index in mpc.bus
            vmap.tab.status     = 5;    %% ones
            vmap.tab.p          = 3;    %% zeros
            vmap.tab.q          = 3;    %% zeros
        end

        function vals = table_var_values(obj, var_names, mpc)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% initialize variable map (all idx, tab = 0)
            vmap = obj.table_var_map(var_names);

            vals = cell(size(var_names));
            r = find(mpc.bus(:, GS) | mpc.bus(:, BS));
            obj.bus = r;
            nr = size(r, 1);
            for k = 1:length(var_names)
                vn = var_names{k};
                switch vmap.tab.(vn)
                    case 0      %% default 'bus' table
                        vals{k} = mpc.bus(r, vmap.idx.(vn));
                    case 1      %% empty char
                        vals{k} = cell(nr, 1);
                        [vals{k}{:}] = deal('');
                    case 2
                        vals{k} = [1:nr]';
                    case 3
                        vals{k} = zeros(nr, 1);
                    case 4
                        vals{k} = r;
                    case 5
                        vals{k} = ones(nr, 1);
                end
            end
        end
    end     %% methods
end         %% classdef
