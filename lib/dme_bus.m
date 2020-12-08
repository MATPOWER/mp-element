classdef dme_bus < dm_element
%DME_BUS  Abstract base class for MATPOWER data model bus table

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        isref   %% ref bus indicator vector for buses that are on
        ispv    %% PV bus indicator vector for buses that are on
        ispq    %% PQ bus indicator vector for buses that are on
        Vm0     %% initial voltage magnitudes (p.u.) for buses that are on
        Va0     %% initial voltage angles (radians) for buses that are on
        Vmin    %% voltage magnitude lower bounds for buses that are on
        Vmax    %% voltage magnitude upper bounds for buses that are on
    end     %% properties

    methods
        %% constructor
        function obj = dme_bus()
            obj@dm_element();   %% call parent constructor
            obj.name = 'bus';
            obj.table = 'bus';
        end
    end     %% methods
end         %% classdef
