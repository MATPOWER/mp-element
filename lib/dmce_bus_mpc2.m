classdef dmce_bus_mpc2 < dmc_element % & dmce_bus
%DMCE_BUS_MPC2  Data model converter for bus elements for MATPOWER case v2.

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
            name = 'bus';
        end

        function df = data_field(obj)
            df = 'bus';
        end

        function vmap = table_var_map(obj, dme, mpc)
            vmap = table_var_map@dmc_element(obj, dme, mpc);

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
               VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            bni_fcn = @(ob, vn, nr, r, mpc)bus_name_import(ob, vn, nr, r, mpc, 1);
            bne_fcn = @(ob, vn, nr, r, mpc, dme)bus_name_export(ob, vn, nr, r, mpc, dme, 1);
            bsi_fcn = @(ob, vn, nr, r, mpc)bus_status_import(ob, vn, nr, r, mpc, BUS_TYPE);

            %% mapping for each name, default is {'col', []}
            vmap.uid{2}         = BUS_I;
            vmap.name           = {'fcn', bni_fcn, bne_fcn};    %% fcns to import/export from/to mpc.bus_name
            vmap.status         = {'fcn', bsi_fcn}; %% fcn w/logic for mpc.bus_types
            vmap.source_uid     = {'cell', ''};     %% empty char
            vmap.base_kv{2}     = BASE_KV;
            vmap.type{2}        = BUS_TYPE;
            vmap.area{2}        = BUS_AREA;
            vmap.zone{2}        = ZONE;
            vmap.vm_lb{2}       = VMIN;
            vmap.vm_ub{2}       = VMAX;
            vmap.va{2}          = VA;
            vmap.vm{2}          = VM;
            if isfield(vmap, 'lam_p')
                vmap.lam_p{2}       = LAM_P;
                vmap.lam_q{2}       = LAM_Q;
                vmap.mu_vm_lb{2}    = MU_VMIN;
                vmap.mu_vm_ub{2}    = MU_VMAX;
            end
        end

        function vals = bus_name_import(obj, vn, nr, r, mpc, c)
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

        function mpc = bus_name_export(obj, vn, nr, r, mpc, dme, c)
            bus_names = dme.tab.name;
            if ~all(cellfun(@isempty, bus_names))
                if nr && isempty(r)
                    mpc.bus_name(:, c) = bus_names;
                else
                    mpc.bus_name(r, c) = bus_names;
                end
            end
        end

        function vals = bus_status_import(obj, vn, nr, r, mpc, c)
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
