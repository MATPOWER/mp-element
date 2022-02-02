classdef mp_math < mp_element_container & opt_model
%MP_MATH  MATPOWER mathematical model abstract base class.
%   ?
%
%   MP_MATH provides properties and methods related to the specific
%   problem specification being solved (e.g. power flow, continuation
%   power flow, optimal power flow, etc.) ...
%
%   Properties
%       ? - ?
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        aux_data    %% struct of auxiliary data relevant to the model,
                    %% e.g. can be passed to model constraint functions
    end

    methods
        function obj = build(obj, nm, dm, mpopt)
            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% init_set_types() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  Initialize MP_IDX_MANAGER fields here, if needed,
            %%              after object construction, but before object use.
            if isempty(obj.var)         %% only if not already initialized
                obj.init_set_types();
            end

            %% create element objects for each class with data
            obj.elements = mp_mapped_array();
            for c = obj.element_classes
                mme = c{1}();       %% element constructor
                if is_index_name(nm.elements, mme.name) %% nm element exists
                    obj.elements.add_elements(mme, mme.name);
                end
            end
        end

        function display(obj)
            fprintf('MATH MODEL CLASS : %s\n', class(obj));
            display@opt_model(obj)

            %% elements
            fprintf('\nELEMENTS\n')
            fprintf('========\n')
            fprintf('  name           class\n');
            fprintf(' -------------  --------------------\n');
            for k = 1:length(obj.elements)
                mme = obj.elements{k};
                fprintf('  %-10s     %s\n', mme.name, class(mme));
            end
        end

        function obj = add_aux_data(obj, nm, dm, mpopt)
            %% create aux_data struct
            obj.aux_data = obj.base_aux_data(nm, dm, mpopt);
        end

        function ad = base_aux_data(obj, nm, dm, mpopt)
            %% get model variables
            vvars = nm.model_vvars();
            zvars = nm.model_zvars();
            vars = {vvars{:} zvars{:}};
            vals = cellfun(@(x)nm.params_var(x), vars, 'UniformOutput', false);

            %% get node types
            [ref, pv, pq, by_elm] = nm.node_types(nm, dm);

            %% create aux_data struct
            ad = cell2struct(...
                {vals{:}, {}, ref, pv, pq, ...
                    length(ref), length(pv), length(pq), by_elm}, ...
                {vars{:}, 'var_map', 'ref', 'pv', 'pq', ...
                    'nref', 'npv', 'npq', 'node_type_by_elm'}, 2);
        end

        function obj = add_vars(obj, nm, dm, mpopt)
            obj.add_system_vars(nm, dm, mpopt);

            %% each element adds its OPF variables
            for k = 1:length(obj.elements)
                obj.elements{k}.add_vars(obj, nm, dm, mpopt);
            end
        end

        function obj = add_system_vars(obj, nm, dm, mpopt)
        end

        function obj = add_constraints(obj, nm, dm, mpopt)
            obj.add_system_constraints(nm, dm, mpopt);

            %% each element adds its OPF variables
            for k = 1:length(obj.elements)
                obj.elements{k}.add_constraints(obj, nm, dm, mpopt);
            end
        end

        function obj = add_system_constraints(obj, nm, dm, mpopt)
            %% node balance constraints
            obj.add_node_balance_constraints(nm, dm, mpopt);
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
        end

        function obj = add_costs(obj, nm, dm, mpopt)
            obj.add_system_costs(nm, dm, mpopt);

            %% each element adds its OPF variables
            for k = 1:length(obj.elements)
                obj.elements{k}.add_costs(obj, nm, dm, mpopt);
            end
        end

        function obj = add_system_costs(obj, nm, dm, mpopt)
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            opt = struct();
        end
    end     %% methods
end         %% classdef
