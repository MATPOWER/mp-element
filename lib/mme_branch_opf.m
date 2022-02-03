classdef mme_branch_opf < mme_branch

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end
    
    methods
        function [mu_vad_lb, mu_vad_ub] = opf_branch_ang_diff_prices(obj, mm, nme)
            %% shadow prices on angle difference limits
            iang = mm.userdata.ang_diff_constrained_branch_idx;
            mu_vad_lb = zeros(nme.nk, 1);
            mu_vad_ub = mu_vad_lb;
            if length(iang)
                ll = mm.get_idx('lin');
                lambda = mm.soln.lambda;
                mu_vad_lb(iang) = lambda.mu_l(ll.i1.ang:ll.iN.ang);
                mu_vad_ub(iang) = lambda.mu_u(ll.i1.ang:ll.iN.ang);
            end
        end
    end     %% methods
end         %% classdef
