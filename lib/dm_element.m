classdef dm_element < handle
%DM_ELEMENT  Abstract base class for MATPOWER data model elements
%   All concrete DM_ELEMENT objects are expected to inherit from a
%   specific DM_FORMAT_* class as well.

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        tab             %% main data table
        nr              %% total number of rows in table
        n               %% number of online elements
        ID              %% nr x 1 vector of unique IDs, maps dmi to ID
        ID2i            %% max(ID) x 1 vector, maps IDs to row indices
        on              %% n x 1 vector of row indices of online elements
        off             %% (nr-n) x 1 vector of row indices of offline elements
        i2on            %% nr x 1 vector mapping row index to index in on/off respectively
    end     %% properties

    methods
        function name = name(obj)
            name = '';      %% e.g. 'bus', 'gen'
        end

        function name = cxn_type(obj)
            %% char array or cell array of char arrays with names of types of
            %% junction elements, i.e. node-creating elements (e.g. 'bus'),
            %% this element can connect to
            name = '';
        end

        function name = cxn_idx_prop(obj)
            %% char array or cell array of char arrays with names of properties
            %% containing indices of junction elements that define connections
            %% (e.g. {'bus_fr', 'bus_to'})
            name = '';
        end

        function name = cxn_type_prop(obj)
            %% char array or cell array of char arrays (dim must match
            %% cxn_idx_prop) with names of properties containing type of
            %% junction elements for each connection, only used if the
            %% junction element type can vary by element (e.g. some lines
            %% connect to one kind of bus, some to another kind)
            name = '';
        end

        function var_names = table_var_names(obj)
            var_names = {'uid', 'name', 'status', 'source_uid'};
        end

        function vars = export_vars(obj, task)
            vars = {};
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
        end

        function nr = count(obj, dm)
            if isempty(obj.tab)
                nr = 0;
            else
                nr = size(obj.tab, 1);
            end
            obj.nr = nr;
        end

        function obj = initialize(obj, dm)
            %% set the IDs
            ID = obj.get_ID(dm);

            %% create ID --> index mapping
            nr = obj.nr;
            maxID = max(ID);
            if ~isempty(ID)
                if nr > 5000 && maxID > 10 * nr %% use sparse map (save memory)
                    ID2i = sparse(ID, ones(nr,1), 1:nr, maxID, 1);
                else                            %% use dense map (faster)
                    ID2i = accumarray([ID ones(nr,1)], 1:nr, [maxID 1]);
                end
            else
                ID2i = [];
            end
            obj.ID2i = ID2i;

            %% initial on/off status
            obj.init_status(dm);
        end

        function ID = get_ID(obj, dm)
            ID = obj.tab.uid;
            obj.ID = ID;
        end

        function obj = init_status(obj, dm)
        end

        function obj = update_status(obj, dm)
            status = obj.tab.status;
            if ~isempty(status)
                obj.on  = find(  status );
                obj.off = find( ~status );
                obj.n   = length(obj.on);
                obj.i2on = zeros(obj.nr, 1);
                obj.i2on(obj.on ) = (1:obj.n);
                obj.i2on(obj.off) = (1:length(obj.off));
            end
        end

        function obj = build_params(obj, dm)
        end

        function obj = rebuild(obj, dm)
            obj.count(dm);
            obj.initialize(dm);
            obj.init_status(dm);
            obj.build_params(dm);
        end

        function display(obj)
%             if have_feature('octave')
%                 struct(obj)
%             else
%                 display@handle(obj)
%             end
            fprintf('DATA MODEL ELEMENT NAME  : %s\n', obj.name);
            fprintf('DATA MODEL ELEMENT CLASS : %s\n', class(obj));
            fprintf('    # OF ROWS            : %d\n', obj.nr);
            fprintf('    # OF ONLINE ELEMENTS : %d\n', obj.n);
            disp(obj.tab);
        end
    end     %% methods
end         %% classdef
