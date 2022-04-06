classdef dmc_element < handle
%DMC_ELEMENT  Abstract base class for data model converter for indv elements

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function name = name(obj)
            name = '';      %% e.g. 'bus', 'gen'
        end

        function dme = data_model_element(obj, dm, name)
            if nargin < 3
                name = obj.name;
            end
            dme = dm.elements.(name);
        end

        function df = data_field(obj)
            df = '';        %% name of field in d for default data table
        end

        function s = data_subs(obj)
            s = struct('type', '.', 'subs', obj.data_field());
        end

        function TorF = data_exists(obj, d)
            TorF = isfield(d, obj.data_field());
        end

        function [nr, nc, r] = get_import_size(obj, d)
            if obj.data_exists(d)
                %% use size of default table
                [nr, nc] = size(subsref(d, obj.data_subs()));
            else
                [nr, nc] = deal(0);
            end
            r = [];                         %% all rows
        end

        function [nr, nc, r] = get_export_size(obj, dme)
            [nr, nc] = size(dme.tab);   %% use size of default table
            r = [];                     %% all rows
        end

        function vmap = table_var_map(obj, dme, d)
            %% initialize with vmap.(<name>) = {'col', []}, for all <name>
            names = dme.table_var_names();
            vals = cell(size(names));
            [vals{:}] = deal({'col', []});
            vmap = cell2struct(vals, names, 2);
        end

        function dme = import(obj, dme, d)
            %% main table
            names = dme.main_table_var_names();
            vmap = obj.table_var_map(dme, d);
            vals = obj.import_table_values(names, vmap, d);
            if ~isempty(vals)
                table_class = mp_table_class();
                dme.tab = table_class(vals{:}, 'VariableNames', names);
            end
        end

        function vals = import_table_values(obj, var_names, vmap, d)
            ss = obj.data_subs();       %% subscripts for data in d
            nv = length(var_names);     %% number of variables

            %% rows in data table and cols in subsref(d, obj.data_subs)
            [nr, nc, r] = obj.get_import_size(d);

            if nr
                %% initialize variable values
                vals = cell(size(var_names));

                %% assign variable values
                for k = 1:nv
                    vn = var_names{k};
                    vm = vmap.(vn);
                    switch vm{1}    %% switch on type of mapping
                        case {'col'}     %% column of default table
                            vals{k} = obj.import_col(d, ss, vn, nr, nc, r, vm{2:end});
                        case 'num'      %% scalar numerical value
                            vals{k} = vm{2} * ones(nr, 1);
                        case {'cell'}   %% cell array of values
                            vals{k} = cell(nr, 1);
                            [vals{k}{:}] = deal(vm{2});
                        case 'IDs'      %% [1:nr]', consecutive IDs
                            vals{k} = [1:nr]';
                        case 'r'        %% r
                            vals{k} = r;
                        case 'fcn'      %% general function
                            import_fcn = vm{2};
                            vals{k} = import_fcn(obj, vn, nr, r, d);
                        otherwise
                            error('dmc_element/import_table_values: %d is an unknown var map type', vm{1});
                    end
                end

                %% check for unique uid's if not generated
                if ~strcmp(vmap.uid{1}, 'IDs') && strcmp(var_names{1}, 'uid')
                    if length(unique(vals{1})) ~= nr
                        error('dmc_element/import_table_values: ''uid'' values must be unique\ndata contains only %d unique ''uid'' value(s) for %d ''%s'' elements\n', ...
                            length(unique(vals{1})), nr, obj.name);
                    end
                end
            else            %% table does not exist or is empty
                vals = [];
            end
        end

        function vals = import_col(obj, d, ss, vn, nr, nc, r, c, sf)
            if nargin > 8
                if isa(sf, 'function_handle')
                    sf = sf(obj, vn);
                end
            else
                sf = 1;
            end

            %% default to zeros if mpc table does not have this many columns
            if c > nc
                vals = zeros(nr, 1);
            else
                data = subsref(d, ss);
                if nr && isempty(r)
                    vals = sf * data(:, c);
                else
                    vals = sf * data(r, c);
                end
            end
        end

        function d = export(obj, dme, d, var_names, idx)
            if nargin < 4
                var_names = 'all';
            end
            if ischar(var_names)
                if strcmp(var_names, 'all')
                    var_names = dme.table_var_names();
                else
                    var_names = {var_names};
                end
            end

            ss = obj.data_subs();       %% subscripts for data in d
            nv = length(var_names);     %% number of variables
            if nv
                %% rows in data table, cols in subsref(d, obj.data_subs)
                [nr, nc, r] = obj.get_export_size(dme);

                %% get variable map
                vmap = obj.table_var_map(dme, d);
            end

            for k = 1:nv
                vn = var_names{k};
                vm = vmap.(vn);
                switch vm{1}    %% switch on type of mapping
                    case {'col'}     %% column of default table
                        d = obj.export_col(dme, d, ss, vn, nr, nc, r, vm{2:end});
                    case {'num', 'cell', 'IDs', 'r'}
                        %% do nothing
                    case 'fcn'      %% general function
                        export_fcn = vm{2};
                        d = export_fcn(obj, vn, nr, r, d, dme);
                    otherwise
                        error('dmc_element/import_table_values: %d is an unknown var map type', vm{1});
                end
            end
        end

        function d = export_col(obj, dme, d, ss, vn, nr, nc, r, c, sf)
            if nargin > 9
                if isa(sf, 'function_handle')
                    sf = sf(obj, vn);
                end
            else
                sf = 1;
            end

            %% default to zeros if mpc table does not have this many columns
            if sf
                ss(end+1).type = '()';
                if nr && isempty(r)
                    ss(end).subs = {':', c};
                else
                    ss(end).subs = {r, c};
                end
                d = subsasgn(d, ss, dme.tab.(vn) / sf);
            end
        end
    end     %% methods
end         %% classdef
