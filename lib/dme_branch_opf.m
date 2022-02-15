classdef dme_branch_opf < dme_branch
%DME_BRANCH_OPF  MATPOWER data model class for branch data

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
            var_names = horzcat( table_var_names@dme_branch(obj), ...
                {'mu_flow_fr_ub', 'mu_flow_to_ub', ...
                 'mu_vad_lb', 'mu_vad_ub'} );
%                 'mu_sm_fr_ub', 'mu_sm_to_ub', ...
%                 'mu_pl_fr_ub', 'mu_pl_to_ub', ...
%                 'mu_cm_fr_ub', 'mu_cm_to_ub', ...
        end

        function vars = export_vars(obj, task)
            vars = horzcat( export_vars@dme_branch(obj), ...
                {'mu_flow_fr_ub', 'mu_flow_to_ub', ...
                 'mu_vad_lb', 'mu_vad_ub'} );
        end
    end     %% methods
end         %% classdef
