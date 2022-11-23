classdef math_model_cpf_accs < mp.math_model_cpf_acc & mp.mm_shared_pfcpf_accs
%MP.MATH_MODEL_CPF_ACCS  MATPOWER mathematical model for continuation power flow (CPF) problem.
%   ?
%
%   MP.MATH_MODEL_CPF_ACCS ... power flow ...
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
        function tag = form_tag(obj)
            tag = 'accs';
        end

        function name = form_name(obj)
            name = 'AC-cartesian-power';
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            %% power balance constraints
            ad = obj.aux_data;
            fcn = @(x)node_balance_equations_cpf(obj, x, nm);
            obj.add_nln_constraint({'Pmis', 'Qmis', 'Vmis'}, [ad.npv+ad.npq;ad.npq;ad.npv], 1, fcn, []);
        end
    end     %% methods
end         %% classdef
