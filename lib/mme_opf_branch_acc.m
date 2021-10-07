classdef mme_opf_branch_acc < mme_opf_branch_ac

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
            nme = obj.network_model_element(nm);

            %% call parent
            add_constraints@mme_opf_branch_ac(obj, mm, nm, dm, mpopt);

            %% branch angle difference limits
            [Aang, lang, uang, iang] = ...
                dm.elements.branch.opf_branch_ang_diff_params(...
                    dm, mpopt.opf.ignore_angle_lim);
            nang = length(iang);
            if nang
                fcn_ang = @(xx)ang_diff_fcn(nme, xx, Aang, lang, uang);
                hess_ang = @(xx, lam)ang_diff_hess(nme, xx, lam, Aang);
                mm.add_nln_constraint({'angL', 'angU'}, [nang;nang], 0, fcn_ang, hess_ang, {'Vr', 'Vi'});
            end
            mm.userdata.ang_diff_constrained_branch_idx = iang;
        end
    end     %% methods
end         %% classdef
