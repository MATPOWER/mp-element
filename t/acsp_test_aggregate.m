classdef acsp_test_aggregate < acsp_aggregate

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        %% constructor
        function obj = acsp_test_aggregate(varargin)
            obj@acsp_aggregate(varargin{:});
            obj.element_classes{end+1} = @ac_gizmo;
            if isempty(obj.node)    %% skip if constructed from existing object
                obj.init_set_types();   %% should be called in mp_idx_manager
                                        %% constructor, if not for:
                                        %% https://savannah.gnu.org/bugs/?52614
            end
        end
    end     %% methods
end         %% classdef
