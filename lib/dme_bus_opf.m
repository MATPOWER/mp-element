classdef dme_bus_opf < dme_bus
%DME_BUS_OPF  MATPOWER data model class for bus data

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
    end     %% properties

    methods
        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dme_bus(obj), ...
                {'lam_p', 'lam_q', 'mu_vm_lb', 'mu_vm_ub'});
        end

        function vars = export_vars(obj, task)
            vars = horzcat( export_vars@dme_bus(obj), ...
                {'vm_lb', 'vm_ub', 'lam_p', 'lam_q', ...
                    'mu_vm_lb', 'mu_vm_ub'} );
        end
    end     %% methods
end         %% classdef
