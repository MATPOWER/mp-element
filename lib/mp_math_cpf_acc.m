classdef mp_math_cpf_acc < mp_math_cpf
%MP_MATH_CPF_ACC  MATPOWER mathematical model for continuation power flow (CPF) problem.
%   ?
%
%   MP_MATH_CPF_ACC ... power flow ...
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
        function obj = mp_math_cpf_acc()
            obj@mp_math_cpf();
            obj.element_classes = { @mme_bus_pf_acc, @mme_gen_pf_ac, ...
                @mme_branch_pf_ac, @mme_load_cpf, @mme_shunt_cpf };
        end
    end     %% methods
end         %% classdef
