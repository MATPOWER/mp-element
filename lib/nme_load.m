classdef nme_load < nm_element

%   MATPOWER
%   Copyright (c) 2019-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        %% constructor
        function obj = nme_load()
            obj@nm_element();
            obj.np = 1;             %% this is a 1 port element
        end

        function name = name(obj)
            name = 'load';
        end
    end     %% methods
end         %% classdef
