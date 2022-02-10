classdef mp_math_pf_acpi < mp_math_pf & mm_shared_pfcpf_acpi
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
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
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
            obj.element_classes = { @mme_bus_pf_acp, @mme_gen_pf_ac, ...
                @mme_branch_pf_ac, ...
                @mme_bus3p, @mme_gen3p, @mme_line3p, @mme_xfmr3p, ...
                @mme_buslink_pf_acp };
        end

        function tag = form_tag(obj)
            tag = 'acpi';
        end

        function name = form_name(obj)
            name = 'AC-polar-current';
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            %% power balance constraints
            ad = obj.aux_data;
            npvq = ad.npv+ad.npq;
            fcn = @(x)node_balance_equations(obj, x, nm);
            obj.add_nln_constraint({'Irmis', 'Iimis'}, [npvq;npvq], 1, fcn, []);
        end
    end     %% methods
end         %% classdef
