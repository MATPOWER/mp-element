classdef dmce_line3p_mpc2 < dmc_element_mpc2 % & dmce_line3p
%DMCE_LINE3P_MPC2  Data model converter for 3-phase line elements for MATPOWER case v2.

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
            name = 'line3p';
        end

        function table = table(obj)
            table = 'line3p';
        end

        function vmap = table_var_map(obj, var_names, mpc)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names, mpc);

            %% map type for each name (default mapping is -1)
            vmap.name.type          = 2;    %% empty char
            vmap.source_uid.type    = 2;    %% empty char
            vmap.pl1_fr.type        = 0;    %% zeros
            vmap.ql1_fr.type        = 0;    %% zeros
            vmap.pl2_fr.type        = 0;    %% zeros
            vmap.ql2_fr.type        = 0;    %% zeros
            vmap.pl3_fr.type        = 0;    %% zeros
            vmap.ql3_fr.type        = 0;    %% zeros
            vmap.pl1_to.type        = 0;    %% zeros
            vmap.ql1_to.type        = 0;    %% zeros
            vmap.pl2_to.type        = 0;    %% zeros
            vmap.ql2_to.type        = 0;    %% zeros
            vmap.pl3_to.type        = 0;    %% zeros
            vmap.ql3_to.type        = 0;    %% zeros

            %% map arguments for each name
            vmap.uid.args           = 1;
           %vmap.name.args          = [];
            vmap.status.args        = 4;
           %vmap.source_uid.args    = [];
            vmap.bus_fr.args        = 2;
            vmap.bus_to.args        = 3;
            vmap.lc.args            = 5;
            vmap.len.args           = 6;
        end

        function tab = create_line_construction_table(obj, lc);
            id = lc(:, 1);
            r = lc(:, 2:7);
            x = lc(:, 8:13);
            c = lc(:, 14:19);
            table_class = mp_table_class();
            tab = table_class(id, r, x, c);
%             r11 = lc(:, 2);
%             r21 = lc(:, 3);
%             r31 = lc(:, 4);
%             r22 = lc(:, 5);
%             r32 = lc(:, 6);
%             r33 = lc(:, 7);
%             x11 = lc(:, 8);
%             x21 = lc(:, 9);
%             x31 = lc(:, 10);
%             x22 = lc(:, 11);
%             x32 = lc(:, 12);
%             x33 = lc(:, 13);
%             c11 = lc(:, 14);
%             c21 = lc(:, 15);
%             c31 = lc(:, 16);
%             c22 = lc(:, 17);
%             c32 = lc(:, 18);
%             c33 = lc(:, 19);
%             if have_feature('table')
%                 tab = table(id, r11, r21, r31, r22, r32, r33, ...
%                                 x11, x21, x31, x22, x32, x33, ...
%                                 c11, c21, c31, c22, c32, c33);
%             else
%                 tab = mp_table( id, r11, r21, r31, r22, r32, r33, ...
%                                     x11, x21, x31, x22, x32, x33, ...
%                                     c11, c21, c31, c22, c32, c33);
%             end
        end

        function dme = import(obj, dme, mpc)
            %% call parent
            dme = import@dmc_element_mpc2(obj, dme, mpc);

            if ~isempty(dme.tab)
                %% system frequency
                dme.freq = mpc.freq;

                %% import line construction table
                dme.lc_tab = create_line_construction_table(obj, mpc.lc);
            end
        end
    end     %% methods
end         %% classdef
