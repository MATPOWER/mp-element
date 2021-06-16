classdef dmce_branch_mpc2 < dmc_element_mpc2 % & dmce_branch
%DMCE_BRANCH_MPC2  Data model converter for branch elements for MATPOWER case v2.

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
        function obj = dmce_branch_mpc2()
            obj.name = 'branch';
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
            vmap.idx.status     = BR_STATUS;
           %vmap.idx.source_uid = 0;
            vmap.idx.bus_fr     = F_BUS;
            vmap.idx.bus_to     = T_BUS;
            vmap.idx.r          = BR_R;
            vmap.idx.x          = BR_X;
            vmap.idx.g_fr       = 0;
            vmap.idx.b_fr       = BR_B;
            vmap.idx.g_to       = 0;
            vmap.idx.b_to       = BR_B;
            vmap.idx.sm_ub_a    = RATE_A;
            vmap.idx.sm_ub_b    = RATE_B;
            vmap.idx.sm_ub_c    = RATE_C;
            vmap.idx.cm_ub_a    = RATE_A;
            vmap.idx.cm_ub_b    = RATE_B;
            vmap.idx.cm_ub_c    = RATE_C;
            vmap.idx.vad_lb     = ANGMIN;
            vmap.idx.vad_ub     = ANGMAX;
            vmap.idx.tm         = TAP;
            vmap.idx.ta         = SHIFT;
            vmap.idx.pl_fr      = PF;
            vmap.idx.ql_fr      = QF;
            vmap.idx.pl_to      = PT;
            vmap.idx.ql_to      = QT;

            %% map table for each name (0 for default mapping)
            vmap.tab.uid        = 2;    %% consecutive IDs, starting at 1
            vmap.tab.name       = 1;    %% empty char
            vmap.tab.source_uid = 1;    %% empty char
            vmap.tab.g_fr       = 3;    %% zeros
            vmap.tab.b_fr       = 4;    %% div by 2
            vmap.tab.g_to       = 3;    %% zeros
            vmap.tab.b_to       = 4;    %% div by 2
        end

        function vals = table_var_values(obj, var_names, mpc)
            %% get variable map (all idx, tab = 0)
            vmap = obj.table_var_map(var_names);

            vals = cell(size(var_names));
            [nr, nc] = size(mpc.branch);
            for k = 1:length(var_names)
                vn = var_names{k};
                switch vmap.tab.(vn)
                    case 0      %% default 'branch' table
                        c = vmap.idx.(vn);
                        if c > nc
                            vals{k} = zeros(nr, 1);
                        else
                            vals{k} = mpc.branch(:, c);
                        end
                    case 1      %% empty char
                        vals{k} = cell(nr, 1);
                        [vals{k}{:}] = deal('');
                    case 2
                        vals{k} = [1:nr]';
                    case 3
                        vals{k} = zeros(nr, 1);
                    case 4
                        vals{k} = mpc.branch(:, vmap.idx.(vn)) / 2;
                end
            end
        end
    end     %% methods
end         %% classdef
