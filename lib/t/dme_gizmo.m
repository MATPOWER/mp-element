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
    end     %% methods
end         %% classdef
