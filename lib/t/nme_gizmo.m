classdef nme_gizmo < nm_element

%   MATPOWER
%   Copyright (c) 2019-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end
    
    methods
        %% constructor
        function obj = nme_gizmo()
            obj@nm_element();
            obj.np = 3;             %% this is a 3 port element
            obj.nz = 2;             %% each with 2 state variables
        end

        function name = name(obj)
            name = 'gizmo';
        end
    end     %% methods
end         %% classdef
