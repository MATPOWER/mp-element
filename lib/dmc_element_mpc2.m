classdef dmc_element_mpc2 < dmc_element
%DMC_ELEMENT_MPC2  Base class for MPC2 data model converter for indv elements

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        table   %% mpc field for default table, assigned in subclass constructor
    end     %% properties

    methods
        function [nr, nc, r] = get_size(obj, mpc)
            [nr, nc] = size(mpc.(obj.table));   %% use size of default table
            r = [];                             %% all rows
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
            nv = length(var_names);             %% number of variables
            [nr, nc, r] = obj.get_size(mpc);    %% rows in data table and
                                                %% cols in mpc.(obj.table)
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
                        %% args = @fcn where
                        %%  vals = fcn(obj, vn, nr, r, mpc)
                        vals{k} = args(obj, vn, nr, r, mpc);
                    otherwise
                        error('dmc_element_mpc2/table_var_values: %d is an unknown var map type', vm.type);
                end
            end
        end
    end     %% methods
end         %% classdef
