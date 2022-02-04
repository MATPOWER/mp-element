classdef mp_math_opf_legacy < handle
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
        function obj = def_set_types_legacy(obj)
            obj.set_types = struct(...
                    'var',  'VARIABLES', ...
                    'lin',  'LINEAR CONSTRAINTS', ...
                    'nle',  'NONLIN EQ CONSTRAINTS', ...
                    'nli',  'NONLIN INEQ CONSTRAINTS', ...
                    'qdc',  'QUADRATIC COSTS', ...
                    'nlc',  'GEN NONLIN COSTS', ...
                    'cost', 'LEGACY COSTS'  ...
                );
        end

        function obj = init_set_types_legacy(obj)
            %% finish initializing data structures for each type
            es = struct();  %% empty struct
            obj.cost.data = struct( ...
                'N', es, ...
                'H', es, ...
                'Cw', es, ...
                'dd', es, ...
                'rh', es, ...
                'kk', es, ...
                'mm', es, ...
                'vs', es );
            obj.cost.params = [];
        end

        function obj = build_legacy(obj, nm, dm, mpopt)
            if strcmp(obj.form_tag, 'dc') && toggle_softlims(obj.mpc, 'status')
                %% user data required by toggle_softlims
                branch_nme = nm.elements.branch;
                [Bbr, pbr] = branch_nme.get_params(1:branch_nme.nk, {'B', 'p'});
                obj.userdata.Bf = Bbr * branch_nme.C';
                obj.userdata.Pfinj = pbr;
            end

            %% execute userfcn callbacks for 'formulation' stage
            if isfield(obj.mpc, 'userfcn')
                userfcn = obj.mpc.userfcn;
            else
                userfcn = [];
            end
            obj = run_userfcn(userfcn, 'formulation', obj, mpopt);
        end
    end     %% methods
end         %% classdef
