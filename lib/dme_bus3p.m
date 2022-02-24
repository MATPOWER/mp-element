classdef dme_bus3p < dm_element
%DME_BUS3P  MATPOWER data model class for 3-phase bus data

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        type        %% node type vector for buses that are on
        vm1_start   %% initial phase 1 voltage magnitudes (p.u.) for buses that are on
        vm2_start   %% initial phase 2 voltage magnitudes (p.u.) for buses that are on
        vm3_start   %% initial phase 3 voltage magnitudes (p.u.) for buses that are on
        va1_start   %% initial phase 1 voltage angles (radians) for buses that are on
        va2_start   %% initial phase 2 voltage angles (radians) for buses that are on
        va3_start   %% initial phase 3 voltage angles (radians) for buses that are on
    end     %% properties

    methods
        function name = name(obj)
            name = 'bus3p';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'type', 'base_kv', 'vm1', 'vm2', 'vm3', 'va1', 'va2', 'va3'});
        end

%         function vars = export_vars(obj, task)
%             switch task
%                 case 'PF'
%                     vars = {'type', 'vm', 'va'};
%                 case 'CPF'
%                     vars = {'type', 'vm', 'va'};
%                 case 'OPF'
%                     vars = {'vm', 'va', 'vm_lb', 'vm_ub', 'lam_p', 'lam_q', ...
%                         'mu_vm_lb', 'mu_vm_ub'};
%                 otherwise
%                     vars = 'all';
%             end
%         end

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

        function vm_start = set_vm_start(obj, p, gen_dme, gbus, ig)
            gen = gen_dme.tab;
            vm_start = obj.tab.(sprintf('vm%d', p))(obj.on);

            %% pull PV bus voltage magnitudes from gen.vm_setpoint
            vcb = ones(obj.n, 1);   %% create mask of voltage-controlled buses
            vcb(obj.type == NODE_TYPE.PQ) = 0;  %% exclude PQ buses
            %% find indices of online gens at online v-c buses
            k = find(vcb(gbus));
            vm_setpoint_prop = sprintf('vm%d_setpoint', p);
            if isempty(ig)
                vm_start(gbus(k)) = gen.(vm_setpoint_prop)(gen_dme.on(k));
            else
                vm_start(gbus(k)) = gen.(vm_setpoint_prop)(gen_dme.on(ig(k)));
            end
        end

        function obj = build_params(obj, dm)
            %% initialize voltage from bus table
            bus = obj.tab;
            obj.va1_start = bus.va1(obj.on) * pi/180;
            obj.va2_start = bus.va2(obj.on) * pi/180;
            obj.va3_start = bus.va3(obj.on) * pi/180;
            obj.vm1_start = bus.vm1(obj.on);
            obj.vm2_start = bus.vm2(obj.on);
            obj.vm3_start = bus.vm3(obj.on);

            %% update bus type and starting vm based on connected gen status
            if dm.elements.is_index_name('gen3p')
                gen_dme = dm.elements.gen3p;
                [gbus, ig] = obj.gbus_vector(gen_dme);
                nb = obj.n;
                ng = length(gbus);

                %% update bus types based on connected generator status
                %% gen connection matrix, element i, j is 1 if gen j @ bus i is ON
                Cg = sparse(gbus, (1:ng)', 1, nb, ng);
                bus_gen_status = Cg * ones(ng, 1);  %% num of gens ON at each bus
%                 obj.type(obj.type == NODE_TYPE.REF & ~bus_gen_status) = NODE_TYPE.PQ;
                  % above line would affect OPF (not just PF, CPF) where REF is
                  % used only as angle reference and does not require an online gen
                obj.type(obj.type == NODE_TYPE.PV & ~bus_gen_status) = NODE_TYPE.PQ;
%                 obj.ensure_ref_bus();   %% pick a new ref bus if one does not exist

                obj.vm1_start = obj.set_vm_start(1, gen_dme, gbus, ig);
                obj.vm2_start = obj.set_vm_start(2, gen_dme, gbus, ig);
                obj.vm3_start = obj.set_vm_start(3, gen_dme, gbus, ig);
            end
        end

%         function obj = set_bus_type_ref(obj, dm, idx)
%             obj.tab.type(obj.on(idx)) = NODE_TYPE.REF;
%             obj.type(idx) = NODE_TYPE.REF;
%         end
% 
%         function obj = set_bus_type_pv(obj, dm, idx)
%             obj.tab.type(obj.on(idx)) = NODE_TYPE.PV;
%             obj.type(idx) = NODE_TYPE.PV;
%         end
% 
%         function obj = set_bus_type_pq(obj, dm, idx)
%             obj.tab.type(obj.on(idx)) = NODE_TYPE.PQ;
%             obj.type(idx) = NODE_TYPE.PQ;
%         end
    end     %% methods
end         %% classdef
