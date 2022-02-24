classdef dmc_element_mpc2 < dmc_element
%DMC_ELEMENT_MPC2  Base class for MPC2 data model converter for indv elements

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
        function table = table(obj)
            table = '';     %% mpc field for default table
        end

        function [nr, nc, r] = get_import_size(obj, mpc)
            if isfield(mpc, obj.table)
                [nr, nc] = size(mpc.(obj.table));   %% use size of default table
            else
                nr = 0; nc = 0;
            end
            r = [];                             %% all rows
        end

        function [nr, nc, r] = get_export_size(obj, dme)
            [nr, nc] = size(dme.tab);   %% use size of default table
            r = [];                     %% all rows
        end

        function vmap = table_var_map(obj, var_names, mpc)
            %% vmap.(name).type = -1, 0, 1, 2, 3, 4, 5 (-1 is default)
            %% vmap.(name).args = <anything>, depends on type

            %% initialize with vmap.(<name>).type = -1
            %%             and vmap.(<name>).args = [], for all <name>
            vals = cell(size(var_names));
            [vals{:}] = deal(struct('type', -1, 'args', []));
            vmap = cell2struct(vals, var_names, 2);
        end

        function vals = table_var_values(obj, var_names, mpc)
            nv = length(var_names);                 %% number of variables
            [nr, nc, r] = obj.get_import_size(mpc); %% rows in data table and
                                                    %% cols in mpc.(obj.table)

            if nr
                %% get variable map
                vmap = obj.table_var_map(var_names, mpc);

                %% initialize variable values
                vals = cell(size(var_names));

                %% assign variable values
                for k = 1:nv
                    vn = var_names{k};
                    vm = vmap.(vn);
                    switch vm.type      %% switch on type of mapping
                        case -1     %% column of default table
                            %% vm.args: c - column index
                            %%          {c, sf} - col idx, scale factor
                            %%              sf = constant or @fcn where
                            %%                  sf = fcn(obj, vn)
                            if iscell(vm.args)      %% get column index
                                c = vm.args{1};
                                if isa(vm.args{2}, 'function_handle')
                                    sf = vm.args{2}(obj, vn);   %% scalar constant
                                else
                                    sf = vm.args{2};            %% scalar constant
                                end
                            else
                                c = vm.args;
                                sf = 1;
                            end
                            if c > nc   %% default to zeros if mpc table does not
                                        %% have this many columns
                                vals{k} = zeros(nr, 1);
                            else
                                if nr && isempty(r)
                                    vals{k} = sf * mpc.(obj.table)(:, c);
                                else
                                    vals{k} = sf * mpc.(obj.table)(r, c);
                                end
                            end
                        case 0      %% zeros
                            vals{k} = zeros(nr, 1);
                        case 1      %% ones
                            vals{k} = ones(nr, 1);
                        case 2      %% empty char
                            vals{k} = cell(nr, 1);
                            [vals{k}{:}] = deal('');
                        case 3      %% [1:nr]', consecutive IDs
                            vals{k} = [1:nr]';
                        case 4      %% r
                            vals{k} = r;
                        case 5      %% general function
                            %% vm.args = @import_fcn or {@import_fcn, @export_fcn}
                            %%  where
                            %%      vals = import_fcn(obj, vn, nr, r, mpc)
                            %%      mpc = export_fcn(obj, vn, nr, r, mpc, dme)
                            import_fcn = vm.args;
                            if iscell(import_fcn)
                                import_fcn = import_fcn{1};
                            end
                            vals{k} = import_fcn(obj, vn, nr, r, mpc);
                        otherwise
                            error('dmc_element_mpc2/table_var_values: %d is an unknown var map type', vm.type);
                    end
                end
            else            %% table does not exist or is empty
                vals = [];
            end
        end

        function mpc = export(obj, dme, mpc, var_names, idx)
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

            nv = length(var_names);     %% number of variables
            if nv
                %% rows in data table, cols in mpc.(obj.table)
                [nr, nc, r] = obj.get_export_size(dme); 

                %% get variable map
                vmap = obj.table_var_map(var_names, mpc);
            end

            for k = 1:nv
                vn = var_names{k};
                vm = vmap.(vn);
                switch vm.type
                    case -1     %% column of default table
                        %% vm.args: c - column index
                        %%          {c, sf} - col idx, scale factor
                        %%              sf = constant or @fcn where
                        %%                  sf = fcn(obj, vn)
                        if iscell(vm.args)      %% get column index
                            c = vm.args{1};
                            if isa(vm.args{2}, 'function_handle')
                                sf = vm.args{2}(obj, vn);   %% scalar constant
                            else
                                sf = vm.args{2};            %% scalar constant
                            end
                        else
                            c = vm.args;
                            sf = 1;
                        end
                        if sf
                            if nr && isempty(r)
                                mpc.(obj.table)(:, c) = dme.tab.(vn) / sf;
                            else
                                mpc.(obj.table)(r, c) = dme.tab.(vn) / sf;
                            end
                        end
                    case 5      %% general function
                        if iscell(vm.args) && length(vm.args) > 1
                            export_fcn = vm.args{2};
                            mpc = export_fcn(obj, vn, nr, r, mpc, dme);
                        end
                end
            end
        end
    end     %% methods
end         %% classdef
