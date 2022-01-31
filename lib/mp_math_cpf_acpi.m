classdef mp_math_cpf_acpi < mp_math_cpf_acp & mm_pf_shared_acpi
%MP_MATH_CPF_ACPI  MATPOWER mathematical model for continuation power flow (CPF) problem.
%   ?
%
%   MP_MATH_CPF_ACPI ... power flow ...
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
        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            %% power balance constraints
            ad = obj.aux_data;
            npvq = ad.npv+ad.npq;
            fcn = @(x)cpf_node_balance_equations(obj, x, nm, ad);
            obj.add_nln_constraint({'Irmis', 'Iimis'}, [npvq;npvq], 1, fcn, []);
        end
    end     %% methods
end         %% classdef
