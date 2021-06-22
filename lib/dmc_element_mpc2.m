classdef dmc_element_mpc2 < dmc_element
%DMC_ELEMENT_MPC2  Base class for MPC2 data model converter for indv elements

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
        function vmap = table_var_map(obj, var_names)
            %% vmap.(name).type = 0, 1, 2 ... (2 is default)
            %% vmap.(name).args = <anything>, depends on type

            %% initialize with vmap.(<name>).type = 2
            %%             and vmap.(<name>).args = [], for all <name>
            vals = cell(size(var_names));
            [vals{:}] = deal(struct('type', 0, 'args', 0));
            vmap = cell2struct(vals, var_names, 2);
        end
    end     %% methods
end         %% classdef
