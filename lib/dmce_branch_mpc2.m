classdef dmce_branch_mpc2 < dmc_element_mpc2 % & dmce_branch
%DMCE_BRANCH_MPC2  Data model converter for branch elements for MATPOWER case v2.

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
            name = 'branch';
        end

        function df = data_field(obj)
            df = 'branch';
        end

        function vmap = table_var_map(obj, var_names, mpc)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names, mpc);

            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

            %% map type for each name (default mapping is -1)
            vmap.uid.type           = 3;    %% consecutive IDs, starting at 1
            vmap.name.type          = 2;    %% empty char
            vmap.source_uid.type    = 2;    %% empty char
            vmap.g_fr.type          = 0;    %% zeros
            vmap.g_to.type          = 0;    %% zeros

            %% map arguments for each name
           %vmap.uid.args           = [];
           %vmap.name.args          = [];
            vmap.status.args        = BR_STATUS;
           %vmap.source_uid.args    = [];
            vmap.bus_fr.args        = F_BUS;
            vmap.bus_to.args        = T_BUS;
            vmap.r.args             = BR_R;
            vmap.x.args             = BR_X;
           %vmap.g_fr.args          = [];
            vmap.b_fr.args          = {BR_B, 0.5};
           %vmap.g_to.args          = [];
            vmap.b_to.args          = {BR_B, 0.5};
            vmap.sm_ub_a.args       = RATE_A;
            vmap.sm_ub_b.args       = RATE_B;
            vmap.sm_ub_c.args       = RATE_C;
            vmap.cm_ub_a.args       = RATE_A;
            vmap.cm_ub_b.args       = RATE_B;
            vmap.cm_ub_c.args       = RATE_C;
            vmap.vad_lb.args        = ANGMIN;
            vmap.vad_ub.args        = ANGMAX;
            vmap.tm.args            = TAP;
            vmap.ta.args            = SHIFT;
            vmap.pl_fr.args         = PF;
            vmap.ql_fr.args         = QF;
            vmap.pl_to.args         = PT;
            vmap.ql_to.args         = QT;
            vmap.mu_flow_fr_ub.args = MU_SF;
            vmap.mu_flow_to_ub.args = MU_ST;
            vmap.mu_vad_lb.args     = MU_ANGMIN;
            vmap.mu_vad_ub.args     = MU_ANGMAX;
        end
    end     %% methods
end         %% classdef
