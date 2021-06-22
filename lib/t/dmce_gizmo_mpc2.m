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
           %vmap.uid.args        = 0;
           %vmap.name.args       = 0;
           %vmap.status.args     = 0;
           %vmap.source_uid.args = 0;
            vmap.bus_1.args      = 1;
            vmap.bus_2.args      = 2;
            vmap.bus_3.args      = 3;
            vmap.Y1r.args        = 4;
            vmap.Y1i.args        = 5;
            vmap.Y2r.args        = 6;
            vmap.Y2i.args        = 7;
            vmap.Lr.args         = 8;
            vmap.Li.args         = 9;
            vmap.Ir.args         = 10;
            vmap.Ii.args         = 11;
            vmap.M1r.args        = 12;
            vmap.M1i.args        = 13;
            vmap.M2r.args        = 14;
            vmap.M2i.args        = 15;
            vmap.Nr.args         = 16;
            vmap.Ni.args         = 17;
            vmap.Sr.args         = 18;
            vmap.Si.args         = 19;
            vmap.Zr1.args        = 20;
            vmap.Zi1.args        = 21;
            vmap.Zr2.args        = 22;
            vmap.Zi2.args        = 23;

            %% map table for each name (0 for default mapping)
            vmap.uid.type        = 2;    %% consecutive IDs, starting at 1
            vmap.name.type       = 1;    %% empty char
            vmap.status.type     = 3;    %% ones
            vmap.source_uid.type = 1;    %% empty char
        end

        function vals = table_var_values(obj, var_names, mpc)
            %% get variable map (all idx, tab = 0)
            vmap = obj.table_var_map(var_names);

            vals = cell(size(var_names));
            [nr, nc] = size(mpc.gizmo);
            for k = 1:length(var_names)
                vn = var_names{k};
                switch vmap.(vn).type
                    case 0      %% default 'gizmo' table
                        c = vmap.(vn).args;
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
