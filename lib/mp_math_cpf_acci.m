classdef mp_math_cpf_acci < mp_math_cpf & mm_pf_shared_acci
%MP_MATH_CPF_ACCI  MATPOWER mathematical model for continuation power flow (CPF) problem.
%   ?
%
%   MP_MATH_CPF_ACCI ... power flow ...
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
        function obj = mp_math_cpf_acci()
            obj@mp_math_cpf();
            obj.element_classes = { @mme_bus_pf_acc, @mme_gen_pf_ac, ...
                @mme_branch_pf_ac, @mme_load_cpf, @mme_shunt_cpf, ...
                @mme_bus3p, @mme_gen3p, @mme_line3p, @mme_xfmr3p, ...
                @mme_buslink_pf_acc };
        end

        function tag = form_tag(obj)
            tag = 'acci';
        end

        function name = form_name(obj)
            name = 'AC-cartesian-current';
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            %% power balance constraints
            ad = obj.aux_data;
            npvq = ad.npv+ad.npq;
            fcn = @(x)node_balance_equations_cpf(obj, x, nm);
            obj.add_nln_constraint({'Irmis', 'Iimis', 'Vmis'}, [npvq;npvq;ad.npv], 1, fcn, []);
        end
    end     %% methods
end         %% classdef
