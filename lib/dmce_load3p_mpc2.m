classdef dmce_load3p_mpc2 < dmc_element_mpc2 % & dmce_load3p
%DMCE_LOAD3P_MPC2  Data model converter for 3-phase load elements for MATPOWER case v2.

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus
    end     %% properties

    methods
        function name = name(obj)
            name = 'load3p';
        end

        function table = table(obj)
            table = 'load3p';
        end

        function vmap = table_var_map(obj, var_names, mpc)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names, mpc);

            %% map type for each name (default mapping is -1)
            vmap.name.type          = 2;    %% empty char
            vmap.source_uid.type    = 2;    %% empty char

            %% map arguments for each name
            vmap.uid.args           = 1;
           %vmap.name.args          = [];
            vmap.status.args        = 3;
           %vmap.source_uid.args    = [];
            vmap.bus.args           = 2;
            vmap.pd1.args           = 4;
            vmap.pd2.args           = 5;
            vmap.pd3.args           = 6;
            vmap.pf1.args           = 7;
            vmap.pf2.args           = 8;
            vmap.pf3.args           = 9;
        end
    end     %% methods
end         %% classdef
