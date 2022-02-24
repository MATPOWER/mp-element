classdef nme_branch < nm_element

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
        function obj = nme_branch()
            obj@nm_element();
            obj.np = 2;             %% this is a 2 port element
        end

        function name = name(obj)
            name = 'branch';
        end
    end     %% methods
end         %% classdef
