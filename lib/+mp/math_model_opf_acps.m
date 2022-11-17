classdef math_model_opf_acps < mp.math_model_opf_acp
%MP.MATH_MODEL_OPF_ACPS  MATPOWER mathematical model for AC optimal power flow (OPF) problem.
%   ?
%
%   MP.MATH_MODEL_OPF_ACPS ... power flow ...
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
            tag = 'acps';
        end

        function name = form_name(obj)
            name = 'AC-polar-power';
        end

        function add_node_balance_constraints(obj, nm, dm, mpopt)
            %% power balance constraints
            nn = nm.node.N;             %% number of nodes
            fcn_mis = @(x)obj.nodal_power_balance_fcn(x, nm);
            hess_mis = @(x, lam)obj.nodal_power_balance_hess(x, lam, nm);
            obj.add_nln_constraint({'Pmis', 'Qmis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end

        function [lam_p, lam_q] = node_power_balance_prices(obj, nm)
            %% shadow prices on node power balance
            nne = obj.get_idx('nle');
            lambda = obj.soln.lambda;
            lam_p = lambda.eqnonlin(nne.i1.Pmis:nne.iN.Pmis);
            lam_q = lambda.eqnonlin(nne.i1.Qmis:nne.iN.Qmis);
        end
    end     %% methods
end         %% classdef
