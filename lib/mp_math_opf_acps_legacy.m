classdef mp_math_opf_acps_legacy < mp_math_opf_acps & mp_math_opf_legacy
%MP_MATH_OPF_ACPS_LEGACY  MATPOWER mathematical model for AC optimal power flow (OPF) problem.
%   ?
%
%   MP_MATH_OPF_ACPS_LEGACY ... power flow ...
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
        function obj = mp_math_opf_acps_legacy()
            obj@mp_math_opf_acps();
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
            add_named_set@mp_math_opf_acps(obj, varargin{:});
            obj.add_named_set_legacy(varargin{:});
        end
 
        function obj = def_set_types(obj)
            obj.def_set_types_legacy();
        end

        function obj = init_set_types(obj)
            init_set_types@mp_math_opf_acps(obj);
            obj.init_set_types_legacy();
        end

        function obj = build(obj, nm, dm, mpopt)
            obj.mpc = dm.source;
            build@mp_math_opf_acps(obj, nm, dm, mpopt);
            obj.build_legacy(nm, dm, mpopt);
        end

        function obj = add_vars(obj, nm, dm, mpopt)
            add_vars@mp_math_opf_acps(obj, nm, dm, mpopt);  %% call parent

            %% legacy user-defined variables
            if isfield(dm.userdata, 'legacy_opf_user_mods')
                obj.add_legacy_user_vars(nm, dm, mpopt);
            end
        end

        function add_system_costs(obj, nm, dm, mpopt)
            add_system_costs@mp_math_opf_acps(obj, nm, dm, mpopt);  %% call parent

            %% legacy user-defined costs
            if isfield(dm.userdata, 'legacy_opf_user_mods')
                obj.add_legacy_user_costs(nm, dm, 0);
            end
        end

        function obj = add_system_constraints(obj, nm, dm, mpopt)
            %% call parent
            add_system_constraints@mp_math_opf_acps(obj, nm, dm, mpopt);

            %% legacy user-defined constraints
            if isfield(dm.userdata, 'legacy_opf_user_mods')
                obj.add_legacy_user_constraints_ac(nm, dm, mpopt);
            end
        end

        function names = legacy_user_var_names(obj)
            names = {'Va', 'Vm', 'Pg', 'Qg'};
        end
    end     %% methods
end         %% classdef
