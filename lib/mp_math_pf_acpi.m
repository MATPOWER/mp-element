classdef mp_math_pf_acpi < mp_math_pf & mm_pf_shared_acpi
%MP_MATH_PF_ACPI  MATPOWER mathematical model for AC power flow (PF) problem.
%   ?
%
%   MP_MATH_PF_ACPI ... power flow ...
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
        function obj = mp_math_pf_acpi()
            obj@mp_math_pf();
            obj.element_classes = { @mme_buslink_pf_acp };
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            %% power balance constraints
            ad = obj.aux_data;
            npvq = ad.npv+ad.npq;
            fcn = @(x)pf_node_balance_equations(obj, x, nm, ad);
            obj.add_nln_constraint({'Irmis', 'Iimis'}, [npvq;npvq], 1, fcn, []);
        end
    end     %% methods
end         %% classdef
