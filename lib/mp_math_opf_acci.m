classdef mp_math_opf_acci < mp_math_opf_acc
%MP_MATH_OPF_ACCI  MATPOWER mathematical model for AC optimal power flow (OPF) problem.
%   ?
%
%   MP_MATH_OPF_ACCI ... power flow ...
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
        function add_node_balance_constraints(obj, nm, dm, mpopt)
            %% power balance constraints
            nn = nm.node.N;             %% number of nodes
            fcn_mis = @(x)obj.nodal_current_balance_fcn(x, nm);
            hess_mis = @(x, lam)obj.nodal_current_balance_hess(x, lam, nm);
            obj.add_nln_constraint({'rImis', 'iImis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end
    end     %% methods
end         %% classdef
