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
        end

        function vmap = table_var_map(obj, var_names)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names);

            %% define named indices into data matrices
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            %% map column indices for each name (0 for exceptions)
           %vmap.idx.uid        = 0;
           %vmap.idx.name       = 0;
            vmap.idx.status     = GEN_STATUS;
           %vmap.idx.source_uid = 0;
            vmap.idx.bus        = GEN_BUS;
            vmap.idx.vm_setpoint= VG;
            vmap.idx.pg_lb      = PMIN;
            vmap.idx.pg_ub      = PMAX;
            vmap.idx.qg_lb      = QMIN;
            vmap.idx.qg_ub      = QMAX;
            vmap.idx.pg         = PG;
            vmap.idx.qg         = QG;
            vmap.idx.in_service = GEN_STATUS;

            %% map table for each name (0 for default mapping)
            vmap.tab.uid        = 2;    %% consecutive IDs, starting at 1
            vmap.tab.name       = 1;    %% empty char
            vmap.tab.source_uid = 1;    %% empty char
        end

        function vals = table_var_values(obj, var_names, mpc)
            %% get variable map (all idx, tab = 0)
            vmap = obj.table_var_map(var_names);

            vals = cell(size(var_names));
            nr = size(mpc.gen, 1);
            for k = 1:length(var_names)
                vn = var_names{k};
                switch vmap.tab.(vn)
                    case 0      %% default 'gen' table
                        vals{k} = mpc.gen(:, vmap.idx.(vn));
                    case 1      %% empty char
                        vals{k} = cell(nr, 1);
                        [vals{k}{:}] = deal('');
                    case 2
                        vals{k} = [1:nr]';
                end
            end
        end
    end     %% methods
end         %% classdef
