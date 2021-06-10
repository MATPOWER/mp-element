classdef mp_math_opf_legacy < mp_math_opf
%MP_MATH_OPF_LEGACY  MATPOWER mathematical model for optimal power flow (OPF) problem.
%   ?
%
%   MP_MATH_OPF_LEGACY ... optimal power flow ...
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

    properties
        cost = [];
        mpc = struct();
    end     %% properties

    methods
        %% constructor
        function om = mp_math_opf_legacy(mpc)
            args = {};
            have_mpc = 0;
            if nargin > 0
                if isa(mpc, 'mp_math_opf_legacy')
                    args = { mpc };
                elseif isstruct(mpc)
                    have_mpc = 1;
                end
            end

            %% call parent constructor
            om@mp_math_opf(args{:});

            if have_mpc
                om.mpc = mpc;
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

        function om = def_set_types(om)
            om.set_types = struct(...
                    'var',  'VARIABLES', ...
                    'lin',  'LINEAR CONSTRAINTS', ...
                    'nle',  'NONLIN EQ CONSTRAINTS', ...
                    'nli',  'NONLIN INEQ CONSTRAINTS', ...
                    'qdc',  'QUADRATIC COSTS', ...
                    'nlc',  'GEN NONLIN COSTS', ...
                    'cost', 'LEGACY COSTS'  ...
                );
        end

        function om = init_set_types(om)
            %% call parent to create base data structures for each type
            init_set_types@mp_math_opf(om);

            %% finish initializing data structures for each type
            es = struct();  %% empty struct
            om.cost.data = struct( ...
                'N', es, ...
                'H', es, ...
                'Cw', es, ...
                'dd', es, ...
                'rh', es, ...
                'kk', es, ...
                'mm', es, ...
                'vs', es );
            om.cost.params = [];
        end

        function obj = build(obj, nm, dm, mpopt)
            obj.mpc = dm.mpc;
            build@mp_math_opf(obj, nm, dm, mpopt);

            if strcmp(nm.form_tag, 'dc') && toggle_softlims(obj.mpc, 'status')
                %% user data required by toggle_softlims
                branch_nme = nm.elements.branch;
                [Bbr, pbr] = branch_nme.get_params(1:branch_nme.nk, {'B', 'p'});
                obj.userdata.Bf = Bbr * branch_nme.C';
                obj.userdata.Pfinj = pbr;
            end

            obj = dm.run_userfcn(obj, mpopt);
        end
    end     %% methods
end         %% classdef
