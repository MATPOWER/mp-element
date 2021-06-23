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
            obj.table = 'gizmo';
        end

        function vmap = table_var_map(obj, var_names, mpc)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names);

            %% map type for each name (default mapping is -1)
            vmap.uid.type           = 3;    %% consecutive IDs, starting at 1
            vmap.name.type          = 2;    %% empty char
            vmap.status.type        = 1;    %% ones
            vmap.source_uid.type    = 2;    %% empty char

            %% map arguments for each name
           %vmap.uid.args           = [];
           %vmap.name.args          = [];
           %vmap.status.args        = [];
           %vmap.source_uid.args    = [];
            vmap.bus_1.args         = 1;
            vmap.bus_2.args         = 2;
            vmap.bus_3.args         = 3;
            vmap.Y1r.args           = 4;
            vmap.Y1i.args           = 5;
            vmap.Y2r.args           = 6;
            vmap.Y2i.args           = 7;
            vmap.Lr.args            = 8;
            vmap.Li.args            = 9;
            vmap.Ir.args            = 10;
            vmap.Ii.args            = 11;
            vmap.M1r.args           = 12;
            vmap.M1i.args           = 13;
            vmap.M2r.args           = 14;
            vmap.M2i.args           = 15;
            vmap.Nr.args            = 16;
            vmap.Ni.args            = 17;
            vmap.Sr.args            = 18;
            vmap.Si.args            = 19;
            vmap.Zr1.args           = 20;
            vmap.Zi1.args           = 21;
            vmap.Zr2.args           = 22;
            vmap.Zi2.args           = 23;
        end
    end     %% methods
end         %% classdef
