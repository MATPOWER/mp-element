classdef dc_aggregate < mp_aggregate & dc_model

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end
    
    methods
        %% constructor
        function obj = dc_aggregate(varargin)
%             fprintf('--> dc_aggregate (%d args)\n', nargin);
            obj@mp_aggregate(varargin{:});
            obj.element_classes = { @dc_bus, @dc_gen, @dc_load, @dc_branch };
%             fprintf('<-- dc_aggregate (%d args)\n', nargin);
        end

    end     %% methods
end         %% classdef
