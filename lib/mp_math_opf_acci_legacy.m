classdef mp_math_opf_acci_legacy < mp_math_opf_acci & mp_math_opf_legacy
%MP_MATH_OPF_ACCI_LEGACY  MATPOWER mathematical model for AC optimal power flow (OPF) problem.
%   ?
%
%   MP_MATH_OPF_ACCI_LEGACY ... power flow ...
%
%   Properties
%       ? - ?
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        %% constructor
        function obj = mp_math_opf_acci_legacy()
            obj@mp_math_opf_acci();
            if nargin > 0 && isstruct(mpc)
                obj.mpc = mpc;
            end

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  INIT_SET_TYPES() is called explicitly as needed
            %%              (if om.var is empty) in ADD_VAR(), DISPLAY() and
            %%              INIT_INDEXED_NAME(), after object construction,
            %%              but before object use.
        end

        function obj = add_named_set(obj, varargin)
            % call parent method (also checks for valid type for named set)
            add_named_set@mp_math_opf_acci(obj, varargin{:});
            obj.add_named_set_legacy(varargin{:});
        end
 
        function obj = def_set_types(obj)
            obj.def_set_types_legacy();
        end

        function obj = init_set_types(obj)
            init_set_types@mp_math_opf_acci(obj);
            obj.init_set_types_legacy();
        end

        function obj = build(obj, nm, dm, mpopt)
            obj.mpc = dm.source;
            build@mp_math_opf_acci(obj, nm, dm, mpopt);
            obj.build_legacy(nm, dm, mpopt);
        end
    end     %% methods
end         %% classdef
