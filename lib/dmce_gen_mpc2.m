classdef dmce_gen_mpc2 < dmc_element_mpc2 % & dmce_gen
%DMCE_GEN_MPC2  Data model converter for gen elements for MATPOWER case v2.

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
        function obj = dmce_gen_mpc2()
            obj.name = 'gen';
            obj.table = 'gen';
        end

        function vmap = table_var_map(obj, var_names, mpc)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names);

            %% define named indices into data matrices
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            %% map type for each name (default mapping is -1)
            vmap.uid.type           = 3;    %% consecutive IDs, starting at 1
            vmap.name.type          = 2;    %% empty char
            vmap.source_uid.type    = 2;    %% empty char

            %% map arguments for each name
           %vmap.uid.args           = [];
           %vmap.name.args          = [];
            vmap.status.args        = GEN_STATUS;
           %vmap.source_uid.args    = [];
            vmap.bus.args           = GEN_BUS;
            vmap.vm_setpoint.args   = VG;
            vmap.pg_lb.args         = PMIN;
            vmap.pg_ub.args         = PMAX;
            vmap.qg_lb.args         = QMIN;
            vmap.qg_ub.args         = QMAX;
            vmap.pg.args            = PG;
            vmap.qg.args            = QG;
            vmap.in_service.args    = GEN_STATUS;
        end
    end     %% methods
end         %% classdef
