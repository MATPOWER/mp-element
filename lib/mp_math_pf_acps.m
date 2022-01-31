classdef mp_math_pf_acps < mp_math_pf & mm_pf_shared_acps
%MP_MATH_PF_ACPS  MATPOWER mathematical model for AC power flow (PF) problem.
%   ?
%
%   MP_MATH_PF_ACPS ... power flow ...
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
        function obj = mp_math_pf_acps()
            obj@mp_math_pf();
            obj.element_classes = { @mme_buslink_pf_acp };
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            alg = mpopt.pf.alg;
            ad = obj.aux_data;
            
            %% power balance constraints
            switch alg
                case  {'FDXB', 'FDBX'}
                    fcn = @(x)pf_node_balance_equations(obj, x, nm, ad, 1);
                otherwise
                    fcn = @(x)pf_node_balance_equations(obj, x, nm, ad);
            end
            obj.add_nln_constraint({'Pmis', 'Qmis'}, [ad.npv+ad.npq;ad.npq], 1, fcn, []);
        end
    end     %% methods
end         %% classdef
