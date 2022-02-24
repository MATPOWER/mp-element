classdef dmce_gen3p_mpc2 < dmc_element_mpc2 % & dmce_gen
%DMCE_GEN3P_MPC2  Data model converter for 3-phase gen elements for MATPOWER case v2.

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
            name = 'gen3p';
        end

        function table = table(obj)
            table = 'gen3p';
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
            vmap.vm1_setpoint.args  = 4;
            vmap.vm2_setpoint.args  = 5;
            vmap.vm3_setpoint.args  = 6;
            vmap.pg1.args           = 7;
            vmap.pg2.args           = 8;
            vmap.pg3.args           = 9;
            vmap.pf1.args           = 10;
            vmap.pf2.args           = 11;
            vmap.pf3.args           = 12;
        end
    end     %% methods
end         %% classdef
