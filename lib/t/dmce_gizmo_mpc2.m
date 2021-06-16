classdef dmce_gizmo_mpc2 < dmc_element_mpc2 % & dmce_gizmo
%DMCE_GIZMO_MPC2  Data model converter for gizmo elements for MATPOWER case v2.

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
        function obj = dmce_gizmo_mpc2()
            obj.name = 'gizmo';
        end

        function vmap = table_var_map(obj, var_names)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names);

            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

            %% map column indices for each name (0 for exceptions)
           %vmap.idx.uid        = 0;
           %vmap.idx.name       = 0;
           %vmap.idx.status     = 0;
           %vmap.idx.source_uid = 0;
            vmap.idx.bus_1      = 1;
            vmap.idx.bus_2      = 2;
            vmap.idx.bus_3      = 3;
            vmap.idx.Y1r        = 4;
            vmap.idx.Y1i        = 5;
            vmap.idx.Y2r        = 6;
            vmap.idx.Y2i        = 7;
            vmap.idx.Lr         = 8;
            vmap.idx.Li         = 9;
            vmap.idx.Ir         = 10;
            vmap.idx.Ii         = 11;
            vmap.idx.M1r        = 12;
            vmap.idx.M1i        = 13;
            vmap.idx.M2r        = 14;
            vmap.idx.M2i        = 15;
            vmap.idx.Nr         = 16;
            vmap.idx.Ni         = 17;
            vmap.idx.Sr         = 18;
            vmap.idx.Si         = 19;
            vmap.idx.Zr1        = 20;
            vmap.idx.Zi1        = 21;
            vmap.idx.Zr2        = 22;
            vmap.idx.Zi2        = 23;

            %% map table for each name (0 for default mapping)
            vmap.tab.uid        = 2;    %% consecutive IDs, starting at 1
            vmap.tab.name       = 1;    %% empty char
            vmap.tab.status     = 3;    %% ones
            vmap.tab.source_uid = 1;    %% empty char
        end

        function vals = table_var_values(obj, var_names, mpc)
            %% get variable map (all idx, tab = 0)
            vmap = obj.table_var_map(var_names);

            vals = cell(size(var_names));
            [nr, nc] = size(mpc.gizmo);
            for k = 1:length(var_names)
                vn = var_names{k};
                switch vmap.tab.(vn)
                    case 0      %% default 'gizmo' table
                        c = vmap.idx.(vn);
                        if c > nc
                            vals{k} = zeros(nr, 1);
                        else
                            vals{k} = mpc.gizmo(:, c);
                        end
                    case 1      %% empty char
                        vals{k} = cell(nr, 1);
                        [vals{k}{:}] = deal('');
                    case 2
                        vals{k} = [1:nr]';
                    case 3
                        vals{k} = ones(nr, 1);
                end
            end
        end
    end     %% methods
end         %% classdef
