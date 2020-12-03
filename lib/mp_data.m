classdef mp_data < handle
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
    end     %% methods
end         %% classdef
