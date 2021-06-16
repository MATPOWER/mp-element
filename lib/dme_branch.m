classdef dme_branch < dm_element
%DME_BRANCH  Abstract base class for MATPOWER data model branch table

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        fbus    %% bus index vector for "from" port (port 1) (all branches)
        tbus    %% bus index vector for "to" port (port 2) (all branches)
        R       %% series resistance (p.u.) for branches that are on
        X       %% series reactance (p.u.) for branches that are on
        B       %% series reactance (p.u.) for branches that are on
        tap     %% transformer off-nominal turns ratio for branches that are on
        shift   %% xformer phase-shift angle (radians) for branches that are on
        rate_a  %% long term flow limit (p.u.) for branches that are on
    end     %% properties

    methods
        %% constructor
        function obj = dme_branch()
            obj@dm_element();   %% call parent constructor
            obj.name = 'branch';
            obj.table = 'branch';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus_fr', 'bus_to', 'r', 'x', 'g_fr', 'b_fr', ...
                'g_to', 'b_to', 'sm_ub_a', 'sm_ub_b', 'sm_ub_c', ...
                'cm_ub_a', 'cm_ub_b', 'cm_ub_c', 'vad_lb', 'vad_ub', ...
                'tm', 'ta', ... %% remove these when we separate out xformers
                'pl_fr', 'ql_fr', 'pl_to', 'ql_to'});
        end
    end     %% methods
end         %% classdef
