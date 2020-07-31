classdef acp_branch < ac_branch & mp_model_acp

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function add_opf_constraints(obj, nm, om, mpc, mpopt)
            %% call parent
            add_opf_constraints@ac_branch(obj, nm, om, mpc, mpopt);

            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

            %% branch voltage angle difference limits
            nb = size(mpc.bus, 1);      %% number of buses
            [Aang, lang, uang, iang] = makeAang(mpc.baseMVA, mpc.branch, nb, mpopt);
            om.add_lin_constraint('ang', Aang, lang, uang, {'Va'});
            om.userdata.iang = iang;
        end
    end     %% methods
end         %% classdef
