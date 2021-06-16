classdef dmce_load_mpc2 < dmc_element_mpc2 % & dmce_load
%DMCE_LOAD_MPC2  Data model converter for load elements for MATPOWER case v2.

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus
    end     %% properties

    methods
        function obj = dmce_load_mpc2()
            obj.name = 'load';
        end

        function vmap = table_var_map(obj, var_names)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names);

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% map column indices for each name (0 for exceptions)
           %vmap.idx.uid        = 0;
           %vmap.idx.name       = 0;
           %vmap.idx.status     = 0;
           %vmap.idx.source_uid = 0;
            vmap.idx.bus        = BUS_I;
            vmap.idx.pd         = PD;
            vmap.idx.qd         = QD;
            vmap.idx.pd_i       = PD;
            vmap.idx.qd_i       = QD;
            vmap.idx.pd_y       = PD;
            vmap.idx.qd_y       = QD;
           %vmap.idx.p          = 0;
           %vmap.idx.q          = 0;

            %% map table for each name (0 for default mapping)
            vmap.tab.uid        = 2;    %% consecutive IDs, starting at 1
            vmap.tab.name       = 1;    %% empty char
            vmap.tab.status     = 6;    %% ones
            vmap.tab.source_uid = 5;    %% index in mpc.bus
            vmap.tab.pd         = 3;    %% scaled nominal load
            vmap.tab.qd         = 3;    %% scaled nominal load
            vmap.tab.pd_i       = 3;    %% scaled nominal load
            vmap.tab.qd_i       = 3;    %% scaled nominal load
            vmap.tab.pd_y       = 3;    %% scaled nominal load
            vmap.tab.qd_y       = 3;    %% scaled nominal load
            vmap.tab.p          = 4;    %% zeros
            vmap.tab.q          = 4;    %% zeros
        end

        function vals = table_var_values(obj, var_names, mpc)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% get scale factors for system wide ZIP loads
            sf = sys_wide_zip_loads(obj, mpc);

            %% get variable map (all idx, tab = 0)
            vmap = obj.table_var_map(var_names);

            vals = cell(size(var_names));
            r = find(mpc.bus(:, PD) | mpc.bus(:, QD));
            obj.bus = r;
            nr = size(r, 1);
            for k = 1:length(var_names)
                vn = var_names{k};
                switch vmap.tab.(vn)
                    case 0      %% default 'bus' table
                        vals{k} = mpc.bus(r, vmap.idx.(vn));
                    case 1      %% empty char
                        vals{k} = cell(nr, 1);
                        [vals{k}{:}] = deal('');
                    case 2
                        vals{k} = [1:nr]';
                    case 3      %% scaled nominal load
                        vals{k} = sf.(vn) * mpc.bus(r, vmap.idx.(vn));
                    case 4      %% zeros
                        vals{k} = zeros(nr, 1);
                    case 5
                        vals{k} = num2cell(obj.bus);
                    case 6
                        vals{k} = ones(nr, 1);
                end
            end
        end

        function sf = sys_wide_zip_loads(obj, mpc)
            if isfield(mpc, 'sys_wide_zip_loads')
                pw = mpc.sys_wide_zip_loads.pw;
                qw = mpc.sys_wide_zip_loads.qw;
                if any(size(pw) ~= [1 3])
                    error('dme_load_mpc2/sys_wide_zip_loads: ''exp.sys_wide_zip_loads.pw'' must be a 1 x 3 vector');
                end
                if abs(sum(pw) - 1) > eps
                    error('dme_load_mpc2/sys_wide_zip_loads: elements of ''exp.sys_wide_zip_loads.pw'' must sum to 1');
                end
                if isempty(qw)
                    qw = pw;
                else
                    if any(size(qw) ~= [1 3])
                        error('dme_load_mpc2/sys_wide_zip_loads: ''exp.sys_wide_zip_loads.qw'' must be a 1 x 3 vector');
                    end
                    if abs(sum(qw) - 1) > eps
                        error('dme_load_mpc2/sys_wide_zip_loads: elements of ''exp.sys_wide_zip_loads.qw'' must sum to 1');
                    end
                end
            else
                pw = [1 0 0];
                qw = [1 0 0];
            end
            %% set scale factors
            sf = struct(    'pd',   pw(1),  'qd',   qw(1), ...
                            'pd_i', pw(2),  'qd_i', qw(2), ...
                            'pd_y', pw(3),  'qd_y', qw(3)   );
        end
    end     %% methods
end         %% classdef
