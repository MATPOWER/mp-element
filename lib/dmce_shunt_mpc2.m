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
           %vmap.uid.args        = 0;
           %vmap.name.args       = 0;
           %vmap.status.args     = 0;
           %vmap.source_uid.args = 0;
            vmap.bus.args        = BUS_I;
            vmap.gs.args         = GS;
            vmap.bs.args         = BS;
           %vmap.p.args          = 0;
           %vmap.q.args          = 0;

            %% map table for each name (0 for default mapping)
            vmap.uid.type        = 2;    %% consecutive IDs, starting at 1
            vmap.name.type       = 1;    %% empty char
            vmap.source_uid.type = 4;    %% index in mpc.bus
            vmap.status.type     = 5;    %% ones
            vmap.p.type          = 3;    %% zeros
            vmap.q.type          = 3;    %% zeros
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
                switch vmap.(vn).type
                    case 0      %% default 'bus' table
                        vals{k} = mpc.bus(r, vmap.(vn).args);
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
