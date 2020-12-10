classdef nme_branch_acp < nme_branch_ac & mp_model_acp

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function add_opf_constraints(obj, nm, om, dm, mpopt)
            %% call parent
            add_opf_constraints@nme_branch_ac(obj, nm, om, dm, mpopt);

            %% branch voltage angle difference limits
            [Aang, lang, uang, iang] = ...
                dm.branch_angle_diff_constraint(mpopt.opf.ignore_angle_lim);
            om.add_lin_constraint('ang', Aang, lang, uang, {'Va'});
            om.userdata.iang = iang;
        end
    end     %% methods
end         %% classdef
