classdef mp_network_accs < mp_network_acc & mp_form_accs

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         va = [];
%         vm = [];
%     end
    
    methods
        %%-----  OPF methods  -----
        function [lam_p, lam_q] = opf_node_power_balance_prices(obj, mm)
            %% shadow prices on node power balance
            nne = mm.get_idx('nle');
            lambda = mm.soln.lambda;
            lam_p = lambda.eqnonlin(nne.i1.Pmis:nne.iN.Pmis);
            lam_q = lambda.eqnonlin(nne.i1.Qmis:nne.iN.Qmis);
        end
    end     %% methods
end         %% classdef
