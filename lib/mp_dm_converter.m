classdef mp_dm_converter < mp_element_container
%MP_DM_CONVERTER  Abstract base class for MATPOWER data model converters.

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
        function obj = build(obj)
            %% create element objects for each class
            obj.elements = mp_mapped_array();
            for c = obj.element_classes
                dmce = c{1}();      %% element constructor
                obj.elements.add_elements(dmce, dmce.name);
            end
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

        function dm = import(obj, dm, d)
            for k = 1:length(obj.elements)
                dmce = obj.elements{k};
                dme = dmce.data_model_element(dm);
                dme = dmce.import(dme, d);
            end
        end

        function d = export(obj, dm, d, task)
            if nargin < 4
                task = '';
            end
            for k = 1:length(dm.elements)
                dme = dm.elements{k};
                if obj.elements.is_index_name(dme.name)
                    dmce = obj.elements.(dme.name);
                    vars = dme.export_vars(task);
                    d = dmce.export(dme, d, vars);
                end
            end
        end

        function display(obj)
%             if have_feature('octave')
%                 struct(obj)
%             else
%                 display@handle(obj)
%             end
            fprintf('DATA CONVERTER CLASS : %s\n', class(obj));

            %% elements
            fprintf('\nELEMENTS\n')
            fprintf('========\n')
            fprintf(' i  name          class\n');
            fprintf('-- -----------   --------------------\n');
            for k = 1:length(obj.elements)
                dmce = obj.elements{k};
                fprintf('%2d  %-13s %s\n', k, dmce.name, class(dmce));
%                 dmce
            end
        end
    end     %% methods
end         %% classdef
