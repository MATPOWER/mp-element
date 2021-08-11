classdef dme_bus < dm_element
%DME_BUS  MATPOWER data model class for bus data

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        type    %% node type vector for buses that are on
        Vm0     %% initial voltage magnitudes (p.u.) for buses that are on
        Va0     %% initial voltage angles (radians) for buses that are on
        Vmin    %% voltage magnitude lower bounds for buses that are on
        Vmax    %% voltage magnitude upper bounds for buses that are on
    end     %% properties

    methods
        %% constructor
        function obj = dme_bus()
            obj@dm_element();   %% call parent constructor
            obj.name = 'bus';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'base_kv', 'type', 'area', 'zone', 'vm_lb', 'vm_ub', ...
                'va', 'vm', 'lam_p', 'lam_q', 'mu_vm_lb', 'mu_vm_ub'});
        end

        function vars = export_vars(obj, task)
            switch task
                case 'PF'
                    vars = {'type', 'vm', 'va'};
                case 'CPF'
                    vars = {'type', 'vm', 'va'};
                case 'OPF'
                    vars = {'vm', 'va', 'vm_lb', 'vm_ub', 'lam_p', 'lam_q', ...
                        'mu_vm_lb', 'mu_vm_ub'};
                otherwise
                    vars = 'all';
            end
        end

        function status = get_status(obj, dm)
            %% overrides dm_element/get_status()

            %% check that all buses have a valid type
            bt = obj.tab.type;
            err = find(~NODE_TYPE.is_valid(bt));
            if ~isempty(err)
                error('dme_bus/get_status: bus %d has an invalid type', err);
            end

            %% temporarily set bus type property with dimensions for all buses
            %% (reduced for online buses only in update_status())
            obj.type = bt;
            status = (bt ~= NODE_TYPE.NONE);       %% bus status
            obj.status = status;
        end

        function obj = update_status(obj, dm)
            %% call parent to fill in on/off
            update_status@dm_element(obj, dm);

            %% update bus type property to correspond to online buses only
            obj.type = obj.type(obj.on);
        end

        function [gbus, ig] = gbus_vector(obj, gen_dme)
            %% buses of online gens
            gbus = obj.i2on(gen_dme.bus(gen_dme.on));
            ig = [];
        end

        function Vm0 = set_Vm0(obj, gen_dme, gbus, ig)
            gen = gen_dme.tab;
            Vm0 = obj.tab.vm(obj.on);

            %% pull PV bus voltage magnitudes from gen.vm_setpoint
            vcb = ones(obj.n, 1);   %% create mask of voltage-controlled buses
            vcb(obj.type == NODE_TYPE.PQ) = 0;  %% exclude PQ buses
            %% find indices of online gens at online v-c buses
            k = find(vcb(gbus));
            if isempty(ig)
                Vm0(gbus(k)) = gen.vm_setpoint(gen_dme.on(k));
            else
                Vm0(gbus(k)) = gen.vm_setpoint(gen_dme.on(ig(k)));
            end
        end

        function obj = build_params(obj, dm)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
               MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
               QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            gen_dme = dm.elements.gen;
            [gbus, ig] = obj.gbus_vector(gen_dme);
            nb = obj.n;
            ng = length(gbus);
            bus = obj.tab;

            %% update bus types based on connected generator status
            %% gen connection matrix, element i, j is 1 if gen j @ bus i is ON
            Cg = sparse(gbus, (1:ng)', 1, nb, ng);
            bus_gen_status = Cg * ones(ng, 1);  %% num of gens ON at each bus
%             obj.type(obj.type == NODE_TYPE.REF & ~bus_gen_status) = NODE_TYPE.PQ;
              % above line would affect OPF (not just PF, CPF) where REF is
              % used only as angle reference and does not require an online gen
            obj.type(obj.type == NODE_TYPE.PV & ~bus_gen_status) = NODE_TYPE.PQ;
%             obj.ensure_ref_bus();   %% pick a new ref bus if one does not exist

            %% initialize voltage from bus table
            obj.Va0 = bus.va(obj.on) * pi/180;
            obj.Vm0 = obj.set_Vm0(gen_dme, gbus, ig);
            obj.Vmin = bus.vm_lb(obj.on);
            obj.Vmax = bus.vm_ub(obj.on);
        end

        function obj = set_bus_type_ref(obj, dm, idx)
            obj.tab.type(obj.on(idx)) = NODE_TYPE.REF;
            obj.type(idx) = NODE_TYPE.REF;
        end

        function obj = set_bus_type_pv(obj, dm, idx)
            obj.tab.type(obj.on(idx)) = NODE_TYPE.PV;
            obj.type(idx) = NODE_TYPE.PV;
        end

        function obj = set_bus_type_pq(obj, dm, idx)
            obj.tab.type(obj.on(idx)) = NODE_TYPE.PQ;
            obj.type(idx) = NODE_TYPE.PQ;
        end

        function obj = update(obj, dm, varargin)
            %% obj.update(dm, name1, val1, name2, val2, ...)
            %% obj.update(dm, idx, name1, val1, name2, val2, ...)

            n = length(varargin);
            if rem(n, 2)    %% odd
                idx = obj.on(varargin{1});
                s = 2;      %% starting arg index
            else            %% even
                idx = obj.on;
                s = 1;      %% starting arg index
            end
            for k = s:2:n-1
                val = varargin{k+1};
                switch varargin{k}
                    case 'Va'
                        obj.tab.va(idx) = val * 180/pi;
                    case 'Vm'
                        obj.tab.vm(idx) = val;
                    case 'bus_type'
                        obj.tab.type(idx) = val;
                        obj.type = obj.tab.type(obj.on);
                    case 'lamP'
                        obj.tab.lam_p(idx) = val / dm.baseMVA;
                    case 'lamQ'
                        obj.tab.lam_q(idx) = val / dm.baseMVA;
                    case 'muVmin'
                        obj.tab.mu_vm_lb(idx) = val;
                    case 'muVmax'
                        obj.tab.mu_vm_ub(idx) = val;
                end
            end
        end
    end     %% methods
end         %% classdef
