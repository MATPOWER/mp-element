classdef dmce_bus_ld_mpc2_node_test < dmce_bus_nld_mpc2_node_test % & dmce_bus
%DMCE_BUS_LD_MPC2_NODE_TEST  Data model converter for bus elements for MATPOWER case v2.

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
        function obj = dmce_bus_ld_mpc2_node_test()
            obj.name = 'bus_ld';
        end

        function vmap = table_var_map(obj, var_names)
            vmap = table_var_map@dmce_bus_nld_mpc2_node_test(obj, var_names);

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% map column indices for each name (0 for exceptions)
            vmap.pd.args = PD;
            vmap.qd.args = QD;
            vmap.gs.args = GS;
            vmap.bs.args = BS;

            %% map table for each name (0 for default mapping)
            vmap.source_uid.type = 4;    %% index in mpc.bus
        end

        function vals = table_var_values(obj, var_names, mpc)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% initialize variable map (all idx, tab = 0)
            vmap = obj.table_var_map(var_names);

            vals = cell(size(var_names));
            r = find(mpc.bus(:, PD) | mpc.bus(:, QD) | mpc.bus(:, GS) | mpc.bus(:, BS));
            obj.bus = r;
            nr = size(r, 1);
            for k = 1:length(var_names)
                vn = var_names{k};
                switch vmap.(vn).type
                    case 0      %% default 'bus' table
                        vals{k} = mpc.bus(r, vmap.(vn).args);
                    case 1      %% empty char
                        vals{k} = cell(nr, 1);
                        [vals{k}{:}] = deal('');
                    case 2      %% 'bus_name'
                        if isfield(mpc, 'bus_name')
                            vals{k} = mpc.bus_name(r, vmap.(vn).args);
                        else
                            vals{k} = cell(nr, 1);
                            [vals{k}{:}] = deal('');
                        end
                    case 3      %% status field
                        vals{k} = mpc.bus(r, BUS_TYPE) ~= NONE;
                    case 4
                        vals{k} = num2cell(obj.bus);
                end
            end
        end
    end     %% methods
end         %% classdef
