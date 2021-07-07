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
            obj.table = 'bus';
        end

        function [nr, nc, r] = get_import_size(obj, mpc)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            tab = mpc.(obj.table);
            r = find(tab(:, PD) | tab(:, QD));
            obj.bus = r;
            nr = size(r, 1);
            nc = size(tab, 2);          %% use nc of default table
        end

        function [nr, nc, r] = get_export_size(obj, dme)
            [nr, nc] = size(dme.tab);   %% use size of default table
            r = dme.tab.source_uid;     %% rows in bus matrix
        end

        function vmap = table_var_map(obj, var_names, mpc)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names, mpc);

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% get scale factors for system wide ZIP loads
            zip_sf = obj.sys_wide_zip_loads(mpc);
            sf_fcn = @(ob, vn)scale_factor_fcn(ob, vn, zip_sf);

            %% map type for each name (default mapping is -1)
            vmap.uid.type           = 3;    %% consecutive IDs, starting at 1
            vmap.name.type          = 2;    %% empty char
            vmap.status.type        = 1;    %% ones
            vmap.source_uid.type    = 4;    %% index in mpc.bus
            vmap.p.type             = 0;    %% zeros
            vmap.q.type             = 0;    %% zeros

            %% map arguments for each name
           %vmap.uid.args           = [];
           %vmap.name.args          = [];
           %vmap.status.args        = [];
           %vmap.source_uid.args    = [];
            vmap.bus.args           = BUS_I;
            vmap.pd.args            = {PD, sf_fcn};
            vmap.qd.args            = {QD, sf_fcn};
            vmap.pd_i.args          = {PD, sf_fcn};
            vmap.qd_i.args          = {QD, sf_fcn};
            vmap.pd_z.args          = {PD, sf_fcn};
            vmap.qd_z.args          = {QD, sf_fcn};
           %vmap.p.args             = [];
           %vmap.q.args             = [];
        end

        function sf = scale_factor_fcn(obj, vn, zip_sf)
            sf = zip_sf.(vn);   %% scale factor
        end

        function sf = sys_wide_zip_loads(obj, mpc)
            if isfield(mpc, 'sys_wide_zip_loads')
                pw = mpc.sys_wide_zip_loads.pw;
                qw = mpc.sys_wide_zip_loads.qw;
                if any(size(pw) ~= [1 3])
                    error('dmce_load_mpc2/sys_wide_zip_loads: ''exp.sys_wide_zip_loads.pw'' must be a 1 x 3 vector');
                end
                if abs(sum(pw) - 1) > eps
                    error('dmce_load_mpc2/sys_wide_zip_loads: elements of ''exp.sys_wide_zip_loads.pw'' must sum to 1');
                end
                if isempty(qw)
                    qw = pw;
                else
                    if any(size(qw) ~= [1 3])
                        error('dmce_load_mpc2/sys_wide_zip_loads: ''exp.sys_wide_zip_loads.qw'' must be a 1 x 3 vector');
                    end
                    if abs(sum(qw) - 1) > eps
                        error('dmce_load_mpc2/sys_wide_zip_loads: elements of ''exp.sys_wide_zip_loads.qw'' must sum to 1');
                    end
                end
            else
                pw = [1 0 0];
                qw = [1 0 0];
            end
            %% set scale factors
            sf = struct(    'pd',   pw(1),  'qd',   qw(1), ...
                            'pd_i', pw(2),  'qd_i', qw(2), ...
                            'pd_z', pw(3),  'qd_z', qw(3)   );
        end
    end     %% methods
end         %% classdef
