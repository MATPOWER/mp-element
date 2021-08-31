classdef nme_branch < nm_element

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'branch';
%     end
    
    methods
        %% constructor
        function obj = nme_branch()
            obj@nm_element();
            obj.name = 'branch';
            obj.np = 2;             %% this is a 2 port element
        end

        function [mu_vad_lb, mu_vad_ub] = opf_branch_ang_diff_prices(obj, mm)
            %% shadow prices on angle difference limits
            iang = mm.userdata.ang_diff_constrained_branch_idx;
            mu_vad_lb = zeros(obj.nk, 1);
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
