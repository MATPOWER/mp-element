classdef dmce_bus_mpc2 < dmc_element_mpc2 % & dmce_bus
%DMCE_BUS_MPC2  Data model converter for bus elements for MATPOWER case v2.

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
        function obj = dmce_bus_mpc2()
            obj.name = 'bus';
        end

        function vmap = table_var_map(obj, var_names)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names);

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% map column indices for each name (0 for exceptions)
            vmap.idx.uid        = BUS_I;
            vmap.idx.name       = 1;
            vmap.idx.status     = BUS_TYPE;
           %vmap.idx.source_uid = 0;
            vmap.idx.base_kv    = BASE_KV;
            vmap.idx.type       = BUS_TYPE;
            vmap.idx.area       = BUS_AREA;
            vmap.idx.zone       = ZONE;
            vmap.idx.vm_lb      = VMIN;
            vmap.idx.vm_ub      = VMAX;
            vmap.idx.va         = VA;
            vmap.idx.vm         = VM;

            %% map table for each name (0 for default mapping)
            vmap.tab.name       = 2;    %% mpc.bus_name
            vmap.tab.status     = 3;    %% special logic for mpc.bus_types
            vmap.tab.source_uid = 1;    %% empty char
        end

        function vals = table_var_values(obj, var_names, mpc)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% get variable map (all idx, tab = 0)
            vmap = obj.table_var_map(var_names);

            vals = cell(size(var_names));
            nr = size(mpc.bus, 1);
            for k = 1:length(var_names)
                vn = var_names{k};
                switch vmap.tab.(vn)
                    case 0      %% default 'bus' table
                        vals{k} = mpc.bus(:, vmap.idx.(vn));
                    case 1      %% empty char
                        vals{k} = cell(nr, 1);
                        [vals{k}{:}] = deal('');
                    case 2      %% 'bus_name'
                        if isfield(mpc, 'bus_name')
                            vals{k} = mpc.bus_name(:, vmap.idx.(vn));
                        else
                            vals{k} = cell(nr, 1);
                            [vals{k}{:}] = deal('');
                        end
                    case 3      %% status field
                        vals{k} = mpc.bus(:, BUS_TYPE) ~= NONE;
                end
            end
        end
    end     %% methods
end         %% classdef
