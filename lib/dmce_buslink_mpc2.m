classdef dmce_buslink_mpc2 < dmc_element_mpc2 % & dmce_buslink
%DMCE_BUSLINK_MPC2  Data model converter for 1-to-3-phase buslink elements for MATPOWER case v2.

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
            name = 'buslink';
        end

        function df = data_field(obj)
            df = 'buslink';
        end

        function vmap = table_var_map(obj, var_names, mpc)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names, mpc);

            %% map type for each name (default mapping is -1)
            vmap.name.type          = 2;    %% empty char
            vmap.source_uid.type    = 2;    %% empty char

            %% map arguments for each name
            vmap.uid.args           = 1;
           %vmap.name.args          = [];
            vmap.status.args        = 4;
           %vmap.source_uid.args    = [];
            vmap.bus.args           = 2;
            vmap.bus3p.args         = 3;
        end
    end     %% methods
end         %% classdef
