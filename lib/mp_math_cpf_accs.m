classdef mp_math_cpf_accs < mp_math_cpf & mm_pf_shared_accs
%MP_MATH_CPF_ACCS  MATPOWER mathematical model for continuation power flow (CPF) problem.
%   ?
%
%   MP_MATH_CPF_ACCS ... power flow ...
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
        function obj = mp_math_cpf_accs()
            obj@mp_math_cpf();
            obj.element_classes = { @mme_pf_buslink_acc };
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            %% power balance constraints
            ad = obj.aux_data;
            fcn = @(x)cpf_node_balance_equations(nm, x, ad);
            obj.add_nln_constraint({'Pmis', 'Qmis', 'Vmis'}, [ad.npv+ad.npq;ad.npq;ad.npv], 1, fcn, []);
        end
    end     %% methods
end         %% classdef
