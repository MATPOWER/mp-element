classdef mp_math_cpf < mp_math
%MP_MATH_CPF  MATPOWER math model for continuation power flow (CPF) problem.
%   ?
%
%   MP_MATH_CPF ... continuation power flow ...
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
        function obj = build(obj, nm, dm, mpopt)
            obj.userdata.aux_data = nm.cpf_aux_data(dm, mpopt);
            nm.cpf_add_vars(obj, nm, dm, mpopt);
            nm.cpf_add_constraints(obj, nm, dm, mpopt);
        end
    end     %% methods
end         %% classdef
