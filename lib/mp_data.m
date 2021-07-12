classdef mp_data < mpe_container
%MP_DATA  Abstract base class for MATPOWER data model

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        baseMVA
        userdata = struct();
    end     %% properties

    methods
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
            obj = dmc.import(obj, d);

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

            %% elements
            fprintf('\nELEMENTS\n')
            fprintf('========\n')
            fprintf(' i  name          table        nr         n      class\n');
            fprintf('-- -----------   ---------  --------  --------  --------------------\n');
            for k = 1:length(obj.elements)
                dme = obj.elements{k};
                fprintf('%2d  %-13s %-10s %6d %9d    %s\n', k, dme.name, dme.table, dme.nr, dme.n, class(dme));
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
        function [success, d, pf] = pf_enforce_q_lims(obj, pf, nm, mpopt);
            gen_dme = obj.elements.gen;
            [mn, mx, both] = gen_dme.violated_q_lims(obj, mpopt);

            if ~isempty(both)   %% we have some Q limit violations
                if isempty(mn) && isempty(mx)   %% infeasible
                    if mpopt.verbose
                        fprintf('All %d remaining gens exceed their Q limits : INFEASIBLE PROBLEM\n', length(both));
                    end
                    d = [];
                    success = 0;
                else
                    if mpopt.verbose && ~isempty(mx)
                        fprintf('Gen %d at upper Q limit, converting to PQ bus\n', gen_dme.on(mx));
                    end
                    if mpopt.verbose && ~isempty(mn)
                        fprintf('Gen %d at lower Q limit, converting to PQ bus\n', gen_dme.on(mn));
                    end

                    %% save corresponding limit values
                    pf.fixed_q_qty(mx) = gen_dme.Qmax(mx);
                    pf.fixed_q_qty(mn) = gen_dme.Qmin(mn);
                    mx = [mx;mn];

                    %% set Qg to binding limit
                    gen_dme.update(obj, mx, 'Qg', pf.fixed_q_qty(mx));

                    %% convert to PQ bus
                    bus_dme = obj.elements.bus;
                    ref0 = find(bus_dme.isref);
                    bidx = bus_dme.i2on(gen_dme.bus(gen_dme.on(mx)));   %% bus of mx
                    if length(ref0) > 1 && any(bus_dme.isref(bidx))
                        error('mp_data/pf_enforce_q_lims: Sorry, MATPOWER cannot enforce Q limits for slack buses in systems with multiple slacks.');
                    end
                    %% set bus type to PQ
                    bus_dme.set_bus_type_pq(obj, bidx);
                    %% potentially pick new slack bus
                    ntv = nm.node_types(nm, obj);       %% node type vector
                    [i1, iN] = nm.get_node_idx('bus');  %% bus node indices
                    btv = ntv(i1:iN);                   %% bus type vector

                    %% indicate if there's been a change in slack bus
                    ref = obj.node_type_ref(btv);       %% new ref bus indices
                    if mpopt.verbose && ref ~= ref0
                        fprintf('Bus %d is new slack bus\n', ...
                            bus_dme.ID(bus_dme.on(ref)));
                    end

                    %% save indices to list of Q limited gens
                    pf.fixed_q_idx = [pf.fixed_q_idx; mx];

                    %% set d to the updated mpc for next step
                    obj.mpc = pf.dmc.elements.gen.export( ...
                        obj.elements.gen, obj.mpc, {'pg', 'qg'});
                    d = obj.mpc;
                    success = 1;
                end
            else                %% no more Q violations
                d = [];
                success = 1;
            end
        end
    end     %% methods
end         %% classdef
