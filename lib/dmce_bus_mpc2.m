classdef dmce_bus_mpc2 < dmc_element_mpc2 % & dmce_bus
%DMCE_BUS_MPC2  Data model converter for bus elements for MATPOWER case v2.

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
        function obj = dmce_bus_mpc2()
            obj.name = 'bus';
            obj.table = 'bus';
        end

        function vmap = table_var_map(obj, var_names, mpc)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names);

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% map type for each name (default mapping is -1)
            vmap.name.type          = 5;    %% fcn to import from mpc.bus_name
            vmap.status.type        = 5;    %% fcn w/logic for mpc.bus_types
            vmap.source_uid.type    = 2;    %% empty char

            %% map arguments for each name
            vmap.uid.args           = BUS_I;
            vmap.name.args          = @(ob, vn, nr, r, mpc)import_bus_name(ob, vn, nr, r, mpc, 1);
            vmap.status.args        = @(ob, vn, nr, r, mpc)import_bus_status(ob, vn, nr, r, mpc, BUS_TYPE);
           %vmap.source_uid.args    = [];
            vmap.base_kv.args       = BASE_KV;
            vmap.type.args          = BUS_TYPE;
            vmap.area.args          = BUS_AREA;
            vmap.zone.args          = ZONE;
            vmap.vm_lb.args         = VMIN;
            vmap.vm_ub.args         = VMAX;
            vmap.va.args            = VA;
            vmap.vm.args            = VM;
        end

        function vals = import_bus_name(obj, vn, nr, r, mpc, c)
            if isfield(mpc, 'bus_name')
                if nr && isempty(r)
                    vals = mpc.bus_name(:, c);
                else
                    vals = mpc.bus_name(r, c);
                end
            else
                vals = cell(nr, 1);
                [vals{:}] = deal('');
            end
        end

        function vals = import_bus_status(obj, vn, nr, r, mpc, c)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;

            if nr && isempty(r)
                vals = mpc.bus(:, BUS_TYPE) ~= NONE;
            else
                vals = mpc.bus(r, BUS_TYPE) ~= NONE;
            end
        end
    end     %% methods
end         %% classdef
