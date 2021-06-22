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
            vals = cell(size(var_names));
            [vals{:}] = deal(0);
            vmap.idx = cell2struct(vals, var_names, 2);
            vmap.tab = cell2struct(vals, var_names, 2);
        end
    end     %% methods
end         %% classdef
