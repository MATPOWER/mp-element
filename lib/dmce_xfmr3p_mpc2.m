classdef dmce_xfmr3p_mpc2 < dmc_element_mpc2 % & dmce_branch
%DMCE_XFMR3P_MPC2  Data model converter for 3-phase transformer elements for MATPOWER case v2.

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
        function obj = dmce_xfmr3p_mpc2()
            obj.name = 'xfmr3p';
            obj.table = 'xfmr3p';
        end

        function vmap = table_var_map(obj, var_names, mpc)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names, mpc);

            %% map type for each name (default mapping is -1)
            vmap.name.type          = 2;    %% empty char
            vmap.source_uid.type    = 2;    %% empty char
            vmap.pl1_fr.type        = 0;    %% zeros
            vmap.pf1_fr.type        = 1;    %% ones
            vmap.pl2_fr.type        = 0;    %% zeros
            vmap.pf2_fr.type        = 1;    %% ones
            vmap.pl3_fr.type        = 0;    %% zeros
            vmap.pf3_fr.type        = 1;    %% ones
            vmap.pl1_to.type        = 0;    %% zeros
            vmap.pf1_to.type        = 1;    %% ones
            vmap.pl2_to.type        = 0;    %% zeros
            vmap.pf2_to.type        = 1;    %% ones
            vmap.pl3_to.type        = 0;    %% zeros
            vmap.pf3_to.type        = 1;    %% ones

            %% map arguments for each name
            vmap.uid.args           = 1;
           %vmap.name.args          = [];
            vmap.status.args        = 4;
           %vmap.source_uid.args    = [];
            vmap.bus_fr.args        = 2;
            vmap.bus_to.args        = 3;
            vmap.r.args             = 5;
            vmap.x.args             = 6;
            vmap.base_kva.args      = 7;
            vmap.base_kv.args       = 8;
        end
    end     %% methods
end         %% classdef
