classdef dme_gizmo < dm_element
%DME_GIZMO  Abstract base class for MATPOWER data model gizmo table

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus1        %% bus index vector for port 1
        bus2        %% bus index vector for port 2
        bus3        %% bus index vector for port 3
    end     %% properties

    methods
        %% constructor
        function obj = dme_gizmo()
            obj@dm_element();   %% call parent constructor
            obj.name = 'gizmo';
            obj.table = 'gizmo';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus_1', 'bus_2', 'bus_3', 'Y1r', 'Y1i', 'Y2r', 'Y2i', ...
                'Lr', 'Li', 'Ir', 'Ii', 'M1r', 'M1i', 'M2r', 'M2i', ...
                'Nr', 'Ni', 'Sr', 'Si', 'Zr1', 'Zi1', 'Zr2', 'Zi2'});
        end
    end     %% methods
end         %% classdef
