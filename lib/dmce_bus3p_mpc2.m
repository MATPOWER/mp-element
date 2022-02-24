classdef dmce_bus3p_mpc2 < dmc_element_mpc2 % & dmce_bus
%DMCE_BUS3P_MPC2  Data model converter for 3-phase bus elements for MATPOWER case v2.

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function name = name(obj)
            name = 'bus3p';
        end

        function table = table(obj)
            table = 'bus3p';
        end

        function vmap = table_var_map(obj, var_names, mpc)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names, mpc);

            %% map type for each name (default mapping is -1)
            vmap.status.type        = 5;    %% fcn w/logic for mpc.bus_types
            vmap.name.type          = 2;    %% empty char
            vmap.source_uid.type    = 2;    %% empty char

            bsi_fcn = @(ob, vn, nr, r, mpc)bus_status_import(ob, vn, nr, r, mpc, 2);

            %% map arguments for each name
            vmap.uid.args           = 1;
           %vmap.name.args          = [];
            vmap.status.args        = bsi_fcn;
           %vmap.source_uid.args    = [];
            vmap.base_kv.args       = 3;
            vmap.type.args          = 2;
            vmap.vm1.args           = 4;
            vmap.vm2.args           = 5;
            vmap.vm3.args           = 6;
            vmap.va1.args           = 7;
            vmap.va2.args           = 8;
            vmap.va3.args           = 9;
        end

        function vals = bus_status_import(obj, vn, nr, r, mpc, c)
            if nr
                if isempty(r)
                    vals = mpc.bus3p(:, 2) ~= 4;
                else
                    vals = mpc.bus3p(r, 2) ~= 4;
                end
            else
                vals = [];
            end
        end
    end     %% methods
end         %% classdef
