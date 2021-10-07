classdef mp_data < mp_element_container
%MP_DATA  Base class for MATPOWER data model

%   MATPOWER
%   Copyright (c) 2020-2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        base_mva        %% system per unit MVA base [1]
        base_kva        %% system per unit kVA base [2]
        source          %% source of data (e.g. mpc, MATPOWER case struct)
        userdata = struct();
    end     %% properties
    %% [1]  for balanced single-phase systems/sections, must be provided if
    %%      system includes any 'bus' elements
    %% [2]  for unbalanced 3-phase systems/sections, must be provided if
    %%      system includes any 'bus3p' elements

    methods
        %% constructor
        function obj = mp_data()
            %% call parent constructor
            obj@mp_element_container();
            obj.element_classes = ...
                { @dme_bus, @dme_gen, @dme_load, ...
                    @dme_branch, @dme_shunt, ...
                    @dme_bus3p, @dme_gen3p, @dme_load3p, ...
                    @dme_line3p, @dme_xfmr3p, ...
                    @dme_buslink };
        end

        function new_obj = copy(obj)
            %% make shallow copy of object
            new_obj = eval(class(obj));  %% create new object
            if have_feature('octave')
                s1 = warning('query', 'Octave:classdef-to-struct');
                warning('off', 'Octave:classdef-to-struct');
            end
            props = fieldnames(obj);
            if have_feature('octave')
                warning(s1.state, 'Octave:classdef-to-struct');
            end
            for k = 1:length(props)
                new_obj.(props{k}) = obj.(props{k});
            end

            %% make copies of each individual element
            new_obj.elements = new_obj.elements.copy();
        end

        function obj = build(obj, d, dmc)
            %% create empty element objects for each class
            obj.elements = mp_mapped_array();
            for k = 1:length(obj.element_classes)
                dme_class = obj.element_classes{k};
                dme = dme_class();      %% element constructor
                obj.elements.add_elements(dme, dme.name);
            end

            %% import data from external format
            dmc.import(obj, d);

            %% count element objects for each class
            %% remove if count is zero
            for k = length(obj.elements):-1:1
                obj.elements{k}.count(obj);     %% set nr
                if obj.elements{k}.nr == 0
                    obj.elements.delete_elements(k);
                end
            end

            obj.initialize();       %% get IDs, self status
            obj.update_status();    %% update status wrt connectivity, define on/off
            obj.build_params();     %% get parameters
        end

        function obj = initialize(obj, dm)
            for k = 1:length(obj.elements)
                obj.elements{k}.initialize(obj);
            end
        end

        function obj = update_status(obj, dm)
            for k = 1:length(obj.elements)
                obj.elements{k}.update_status(obj);
            end
        end

        function obj = build_params(obj, dm)
            for k = 1:length(obj.elements)
                obj.elements{k}.build_params(obj);
            end
        end

        function n = online(obj, name)
            if obj.elements.is_index_name(name)
                n = obj.elements.(name).n;
            else
                n = 0;
            end
        end

        function display(obj)
%             if have_feature('octave')
%                 struct(obj)
%             else
%                 display@handle(obj)
%             end
            fprintf('DATA MODEL CLASS : %s\n', class(obj));
            fprintf('        base_mva : %g\n', obj.base_mva);

            %% elements
            fprintf('\nELEMENTS\n')
            fprintf('========\n')
            fprintf(' i  name            nr         n      class\n');
            fprintf('-- -----------   --------  --------  --------------------\n');
            for k = 1:length(obj.elements)
                dme = obj.elements{k};
                fprintf('%2d  %-13s %6d %9d    %s\n', k, dme.name, dme.nr, dme.n, class(dme));
%                 dme
            end

            %% user data
            fields = fieldnames(obj.userdata);
            if ~isempty(fields)
                fprintf('\nUSER DATA\n')
                fprintf('=========\n')
                fprintf('  name                         size       class\n');
                fprintf(' ------------------------   -----------  --------------------\n');
                for k = 1:length(fields)
                    f = obj.userdata.(fields{k});
                    [m, n] = size(f);
                    fprintf('  %-24s  %5dx%-5d   %s\n', fields{k}, m, n, class(f));
                end
            end
        end

        %%-----  PF methods  -----
        function [dm1, dm2] = fdpf_B_matrix_models(obj, alg)
            %% [dmp, dmpp] = obj.fdpf_B_matrix_models(alg)
            %% dmpp = obj.fdpf_B_matrix_models(alg)
            %% returns copies of dm used for building B prime, B double prime
            %% for fast-decoupled power flow

            %% modify data model to form Bp (B prime)
            if nargout > 1      %% for both Bp and Bpp
                dm1 = obj.copy();
                if dm1.elements.is_index_name('shunt')
                    dm1.elements.shunt.tab.bs(:) = 0;   %% zero out shunts at buses
                end
                dm2 = dm1.copy();
                dm1.elements.branch.tab.b_fr(:) = 0;    %% zero out line charging shunts
                dm1.elements.branch.tab.b_to(:) = 0;
                dm1.elements.branch.tab.tm(:) = 1;      %% cancel out taps
                if strcmp(alg, 'FDXB')                  %% if XB method
                    dm1.elements.branch.tab.r(:) = 0;   %% zero out line resistance
                end
                dm1 = dm1.build_params(dm1);
            else
                dm2 = obj.copy();
            end

            %% modify data model to form Bpp (B double prime)
            dm2.elements.branch.tab.ta(:) = 0;      %% zero out phase shifters
            if strcmp(alg, 'FDBX')                  %% if BX method
                dm2.elements.branch.tab.r(:) = 0;   %% zero out line resistance
            end

            if nargout > 1      %% for both Bp and Bpp
                dm2 = dm2.build_params(dm2);
            else                %% for just Bpp
                dm1 = dm2.build_params(dm2);
            end
        end

        %%-----  OPF methods  -----
        function obj = set_bus_v_lims_via_vg(obj, use_vg)
            bus_dme = obj.elements.bus;
            gen_dme = obj.elements.gen;
            gbus = bus_dme.i2on(gen_dme.bus(gen_dme.on));   %% buses of online gens
            nb = bus_dme.n;
            ng = gen_dme.n;

            %% gen connection matrix, element i, j is 1 if, generator j at bus i is ON
            Cg = sparse(gbus, (1:ng)', 1, nb, ng);
            Vbg = Cg * sparse(1:ng, 1:ng, gen_dme.vm_setpoint, ng, ng);
            vm_ub = max(Vbg, [], 2);    %% zero for non-gen buses, else max vm_setpoint of gens @ bus
            ib = find(vm_ub);               %% buses with online gens
            vm_lb = max(2*Cg - Vbg, [], 2); %% same as vm_ub, except min vm_setpoint of gens @ bus
            vm_lb(ib) = 2 - vm_lb(ib);

            if use_vg == 1      %% use vm_setpoint setpoint directly
                bus_dme.vm_ub(ib) = vm_ub(ib);  %% ub set by max vm_setpoint @ bus
                bus_dme.vm_lb(ib) = vm_lb(ib);  %% lb set by min vm_setpoint @ bus
                bus_dme.vm_start(ib) = vm_ub(ib);
            elseif use_vg > 0 && use_vg < 1     %% fractional value
                %% use weighted avg between original vm_lb/vm_ub limits and vm_setpoint
                bus_dme.vm_ub(ib) = (1-use_vg) * bus_dme.vm_ub(ib) + use_vg * vm_ub(ib);
                bus_dme.vm_lb(ib) = (1-use_vg) * bus_dme.vm_lb(ib) + use_vg * vm_lb(ib);
            else
                error('mp_data/set_bus_v_lims_via_vg: option ''opf.use_vg'' (= %g) cannot be negative or greater than 1', use_vg);
            end

            %% update bus table as well (for output)
            bus_dme.tab.vm_ub(bus_dme.on(ib)) = bus_dme.vm_ub(ib);
            bus_dme.tab.vm_lb(bus_dme.on(ib)) = bus_dme.vm_lb(ib);
        end
    end     %% methods
end         %% classdef
