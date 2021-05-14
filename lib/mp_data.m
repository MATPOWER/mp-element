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
    end     %% properties

    methods
        function new_obj = copy(obj)
            %% start with shallow copy of object
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
            for k = 1:length(obj.elm_list)
                new_obj.elm_list{k} = new_obj.elm_list{k}.copy();
            end
        end

        function obj = build(obj, d)
            %% create element objects for each class
            i = 0;
            for c = obj.element_classes
                dme = c{1}();       %% element constructor
                if dme.count(obj)   %% set nr
                    i = i + 1;
                    obj.elm_list{i} = dme;
                    obj.elm_map.(dme.name) = i;
                end
            end

            obj.initialize();       %% get IDs, self status
            obj.update_status();    %% update status wrt connectivity, define on/off
            obj.build_params();     %% get parameters
        end

        function obj = initialize(obj, dm)
            for dme = obj.elm_list
                dme{1}.initialize(obj);
            end
        end

        function obj = update_status(obj, dm)
            for dme = obj.elm_list
                dme{1}.update_status(obj);
            end
        end

        function obj = build_params(obj, dm)
            for dme = obj.elm_list
                dme{1}.build_params(obj);
            end
        end

        function n = online(obj, name)
            dme = obj.elm_by_name(name);
            if isempty(dme)
                n = 0;
            else
                n = dme.n;
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
            fprintf(' i     name        table       nr         n      class\n');
            fprintf('-- -----------   ---------  --------  --------  --------------------\n');
            for k = 1:length(obj.elm_list)
                dme = obj.elm_list{k};
                fprintf('%2d %9s %12s %9d %9d    %s\n', k, dme.name, dme.table, dme.nr, dme.n, class(dme));
%                 dme
            end
        end
    end     %% methods
end         %% classdef
