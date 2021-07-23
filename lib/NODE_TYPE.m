classdef (Sealed) NODE_TYPE
%NODE_TYPE  Defines enumerated type for node types.
%
%   NODE_TYPE.PQ = 1
%   NODE_TYPE.PV = 2
%   NODE_TYPE.REF = 3
%   NODE_TYPE.NONE = 4

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties (Constant)
        PQ = 1;
        PV = 2;
        REF = 3;
        NONE = 4;
    end

%     methods (Access = private)      %% to prevent instantiation
%         function obj = NODE_TYPES
%         end
%     end

    methods (Static)
        function TorF = is_valid(val)
            TorF = val == NODE_TYPE.PQ  | val == NODE_TYPE.PV | ...
                   val == NODE_TYPE.REF | val == NODE_TYPE.NONE;
        end
    end
end
