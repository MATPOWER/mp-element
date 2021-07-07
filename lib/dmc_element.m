classdef dmc_element < handle
%DMC_ELEMENT  Abstract base class for data model converter for indv elements

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        name
    end     %% properties

    methods
        function dme = data_model_element(obj, dm, name)
            if nargin < 3
                name = obj.name;
            end
            dme = dm.elements.(name);
        end

        function dme = import(obj, dme, d)
            var_names = dme.table_var_names();
            var_vals  = obj.table_var_values(var_names, d);
            if have_feature('table')
                dme.tab = table(var_vals{:}, 'VariableNames', var_names);
            else
                dme.tab = mp_table(var_vals{:}, 'VariableNames', var_names);
            end
        end
    end     %% methods
end         %% classdef
