classdef dc_aggregate < mp_aggregate & dc_model

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        v = [];
        z = [];
    end
    
    methods
        %% constructor
        function obj = dc_aggregate(varargin)
%             fprintf('--> dc_aggregate (%d args)\n', nargin);
            obj@mp_aggregate(varargin{:});
            obj.element_classes = { @dc_bus, @dc_gen, @dc_load, @dc_branch };
            if isempty(obj.node)    %% skip if constructed from existing object
                obj.init_set_types();   %% should be called in mp_idx_manager
                                        %% constructor, if not for:
                                        %% https://savannah.gnu.org/bugs/?52614
            end
%             fprintf('<-- dc_aggregate (%d args)\n', nargin);
        end

        function obj = def_set_types(obj)
            def_set_types@mp_aggregate(obj);        %% call parent first
            obj.set_types.v = 'VOLTAGE VARS (v)';
            obj.set_types.z = 'NON-VOLTAGE VARS (z)';
        end

        function obj = build_params(obj, asm, mpc)
            %% call parent to build individual element parameters
            build_params@mp_aggregate(obj, asm, mpc);

            %% aggregate parameters from individual elements
            obj.B = obj.stack_matrix_params('B', 1);
            obj.K = obj.stack_matrix_params('K', 0);
            obj.p = obj.stack_vector_params('p');
        end
    end     %% methods
end         %% classdef
