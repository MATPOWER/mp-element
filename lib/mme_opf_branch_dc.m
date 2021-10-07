classdef mme_opf_branch_dc < mme_branch

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'branch';
%     end
    
    methods
        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %% find branches with flow limits
            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);
            ibr = find(dme.rate_a ~= 0 & dme.rate_a < 1e10);
            nl2 = length(ibr);      %% number of constrained branches
            mm.userdata.flow_constrained_branch_idx = ibr;

            if nl2
                %% limits
                flow_max = dme.rate_a(ibr); %% RATE_A

                %% branch flow constraints
                [B, K, p] = nme.get_params(ibr);
                Af = B * nme.C';
                mm.add_lin_constraint('Pf', Af, -p-flow_max, -p+flow_max, ...
                    nm.va.order);
            end

            %% branch voltage angle difference limits
            [Aang, lang, uang, iang] = ...
                dm.elements.branch.opf_branch_ang_diff_params(...
                    dm, mpopt.opf.ignore_angle_lim);
            if length(iang)
                mm.add_lin_constraint('ang', Aang, lang, uang, {'Va'});
            end
            mm.userdata.ang_diff_constrained_branch_idx = iang;
        end
    end     %% methods
end         %% classdef
