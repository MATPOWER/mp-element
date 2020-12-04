classdef mp_data < handle
%MP_DATA  Abstract base class for MATPOWER data model

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        tab         %% struct array with fields 'name', 'id_col', 'st_col'
                    %% corresponding to the name of the data table, the column
                    %% containing the unique element ID, and the column
                    %% defining the element on/off status, where 0 for id_col
                    %% means the ids are simply the row indices, and 0 for
                    %% st_col means that there is no status column and all
                    %% elements are assumed to be on
        map         %% struct with fields corresponding to OBJ.tab(:).name
                    %% and fields:
                    %%  k      - index into OBJ.tab
                    %%  N      - number of elements
                    %%  n      - number of online elements
                    %%  di2id  - N x 1 vector that maps data model index
                    %%           (e.g. row index in bus table) to unique
                    %%           element ID (e.g. bus number)
                    %%  id2di  - M x 1 vector that maps unique element ID
                    %%           (e.g. bus number) to data model index
                    %%           (e.g. row index in bus table)
                    %%  status - N x 1 on/off status vector
                    %%  on     - n x 1 index vector of online elements
                    %%  off    - (N-n) x 1 index vector of offline elements
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
        end

        function obj = create_mappings(obj)
            for k = 1:length(obj.tab)
                tab = obj.tab(k);

                %% index <--> ID mapping
                [ids, N, M] = obj.get_ids(tab.name, tab.id_col);
                if ~isempty(ids)
                    if N > 5000 && M > 10 * N   %% use sparse map (saves memory)
                        id2di =  sparse(ids, ones(N,1), 1:N, M, 1);
                    else                        %% use dense map (faster)
                        id2di = accumarray([ids ones(N,1)], 1:N, [M 1]);
                    end
                else
                    id2di = [];
                end

                %% initial on/off status
                stats = obj.get_status(tab.name, tab.st_col);

                obj.map.(tab.name) = struct( ...
                        'k',        k, ...      %% index into OBJ.tab
                        'N',        N, ...      %% number of elements
                        'n',        [], ...     %% number of online elements
                        'di2id',    ids, ...    %% vector of unique IDs
                        'id2di',    id2di, ...  %% map from ID to row index
                        'status',   stats, ...  %% on/off status vector
                        'on',       [], ...     %% indices of online elements
                        'off',      [] ...      %% indices of offline elements
                    );
            end
        end

        function obj = update_status(obj)
            %% any final status updates & creation of on/off index vectors
            for k = 1:length(obj.tab)
                tab = obj.tab(k);
                map = obj.map.(tab.name);
                if ~isempty(map.status)
                    on  = find(  map.status );
                    off = find( ~map.status );
                    map.on  = on;
                    map.off = off;
                    map.n   = length(on);
                    obj.map.(tab.name) = map;
                end
            end
        end

        function n = online(obj, tab_name)
            if isfield(obj.map, tab_name)
                n = obj.map.(tab_name).n;
            elseif isfield(obj.mpc, tab_name)
                n = size(obj.mpc.(tab_name), 1);
            else
                n = 0;
            end
        end
    end     %% methods
end         %% classdef
