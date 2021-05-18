classdef mp_math_pf < mp_math
%MP_MATH_PF  MATPOWER mathematical model for power flow (PF) problem.
%   ?
%
%   MP_MATH_PF ... power flow ...
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
            obj.userdata.aux_data = nm.pf_aux_data(dm, mpopt);
            nm.pf_add_vars(obj, nm, dm, mpopt);
            nm.pf_add_constraints(obj, nm, dm, mpopt);
        end
    end     %% methods
end         %% classdef
