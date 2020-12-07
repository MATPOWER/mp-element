classdef dm_element < handle
%DM_ELEMENT  Abstract base class for MATPOWER data model elements
%   All concrete DM_ELEMENT objects are expected to inherit from a
%   specific DM_FORMAT_* class as well.

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        name = 'dm_element';
        table           %% name of table to be checked for presence of this
                        %% element type
        id_col = 0;     %% column containing the unique ID
                        %%      0 means no explicit ID
        st_col = 0;     %% column defining element's on/off status
                        %%      0 means no explicit status (i.e all on)
        nr              %% total number of rows in table
        n               %% number of online elements
        ID              %% nr x 1 vector of unique IDs, maps dmi to ID
                        %% let maxID = max(ID)
        ID2i            %% maxID x 1 vector, maps IDs to row indices
        status          %% nr x 1 vector of on/off status for elements
        on              %% n x 1 vector of row indices of online elements
        off             %% (nr-n) x 1 vector of row indices of offline elements
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

        function obj = create_model(obj, dm)
            %% set the IDs
            ID = obj.get_ID(dm);

            %% create ID --> index mapping
            nr = obj.nr;
            maxID = max(ID);
            if ~isempty(ID)
                if nr > 5000 && maxID > 10 * nr %% use sparse map (save memory)
                    ID2i =  sparse(ID, ones(nr,1), 1:nr, maxID, 1);
                else                            %% use dense map (faster)
                    ID2i = accumarray([ID ones(nr,1)], 1:nr, [maxID 1]);
                end
            else
                ID2i = [];
            end
            obj.ID2i = ID2i;

            %% initial on/off status
            status = obj.get_status(dm);
        end

        function ID = get_ID(obj, dm)
            tab = obj.get_table(dm);
            if obj.id_col
                ID = tab(:, obj.id_col);
            else
                ID = [];
            end
            obj.ID = ID;
        end

        function status = get_status(obj, dm)
            tab = obj.get_table(dm);
            if obj.st_col
                status = tab(:, obj.st_col);
            else
                status = ones(obj.nr, 1);
            end
            obj.status = status;
        end

        function status = update_status(obj, tab)
            if ~isempty(obj.status)
                obj.on  = find(  obj.status );
                obj.off = find( ~obj.status );
                obj.n   = length(obj.on);
            end
        end
    end     %% methods
end         %% classdef
