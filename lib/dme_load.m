classdef dme_load < dm_element
%DME_LOAD  Abstract base class for MATPOWER data model load table

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus     %% bus index vector (all loads)
        Pd      %% active power demand (p.u.) for loads that are on
        Qd      %% reactive power demand (p.u.) for loads that are on
    end     %% properties

    methods
        %% constructor
        function obj = dme_load()
            obj@dm_element();   %% call parent constructor
            obj.name = 'load';
            obj.table = 'bus';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus', 'pd', 'qd', 'pd_i', 'qd_i', 'pd_y', 'qd_y', ...
                'p', 'q'});
        end
    end     %% methods
end         %% classdef
