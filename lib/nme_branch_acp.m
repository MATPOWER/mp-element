classdef nme_branch_acp < nme_branch_ac & mp_form_acp

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function opf_add_constraints(obj, mm, nm, dm, mpopt)
            %% call parent
            opf_add_constraints@nme_branch_ac(obj, mm, nm, dm, mpopt);

            %% branch voltage angle difference limits
            [Aang, lang, uang, iang] = ...
                dm.branch_angle_diff_constraint(mpopt.opf.ignore_angle_lim);
            mm.add_lin_constraint('ang', Aang, lang, uang, {'Va'});
            mm.userdata.iang = iang;
        end
    end     %% methods
end         %% classdef
