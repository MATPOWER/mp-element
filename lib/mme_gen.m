classdef mme_gen < mm_element

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'gen';
%     end
    
    methods
        %% constructor
        function obj = mme_gen()
            obj@mm_element();
            obj.name = 'gen';
        end
    end     %% methods
end         %% classdef
