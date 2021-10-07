classdef mp_math_pf_dc < mp_math_pf & mm_pf_shared_dc
%MP_MATH_PF_DC  MATPOWER mathematical model for DC power flow (PF) problem.
%   ?
%
%   MP_MATH_PF_DC ... power flow ...
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
        function obj = mp_math_pf_dc()
            obj@mp_math_pf();
            obj.element_classes = { };
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            ad = obj.aux_data;
            pvq = [ad.pv; ad.pq];

            %% power balance constraints
            A = ad.B(pvq, pvq);
            b = (ad.Pbus(pvq) - ad.B(pvq, ad.ref) * ad.va(ad.ref));
            obj.add_lin_constraint('Pmis', A, b, b);
        end
    end     %% methods
end         %% classdef
