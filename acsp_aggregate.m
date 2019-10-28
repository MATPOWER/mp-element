classdef acsp_aggregate < mp_aggregate & acsp_model

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
        function obj = acsp_aggregate(varargin)
            obj@mp_aggregate(varargin{:});
            obj.element_classes = ...
                { @acsp_bus, @acsp_gen, @acsp_load, @acsp_branch };
        end

    end     %% methods
end         %% classdef
