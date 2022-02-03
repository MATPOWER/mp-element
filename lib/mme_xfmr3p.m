classdef mme_xfmr3p < mm_element

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'xfmr3p';
%     end
    
    methods
        %% constructor
        function obj = mme_xfmr3p()
            obj@mm_element();
            obj.name = 'xfmr3p';
        end
    end     %% methods
end         %% classdef
