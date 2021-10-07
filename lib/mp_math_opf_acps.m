classdef mp_math_opf_acps < mp_math_opf_ac
%MP_MATH_OPF_ACPS  MATPOWER mathematical model for AC optimal power flow (OPF) problem.
%   ?
%
%   MP_MATH_OPF_ACPS ... power flow ...
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
        function obj = mp_math_opf_acps()
            obj@mp_math_opf_ac();
            obj.element_classes = { @mme_opf_gen_ac, @mme_opf_branch_acp, ...
                @mme_opf_buslink_acp };
        end

        function add_node_balance_constraints(obj, nm, dm, mpopt)
            %% power balance constraints
            nn = nm.node.N;             %% number of nodes
            fcn_mis = @(x)opf_power_balance_fcn(nm, nm.opf_convert_x(x, obj.aux_data));
            hess_mis = @(x, lam)opf_power_balance_hess(nm, ...
                nm.opf_convert_x(x, obj.aux_data), lam);
            obj.add_nln_constraint({'Pmis', 'Qmis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end
    end     %% methods
end         %% classdef
