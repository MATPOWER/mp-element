classdef mpe_branch_acc_nln < mpe_branch_acc & mpe_wrapper_ac_nln

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpe_class = @mpe_branch_acc;
    end

    methods
        function obj = mpe_branch_acc_nln()
            obj@mpe_branch_acc();
            obj.mpe_wrapper_ac_nln_init();
        end

        function build_params(obj, nm, mpc)
            build_params@mpe_branch_acc(obj, nm, mpc);
            obj.build_nln_params(nm, mpc);
        end

        function nk = count(obj, mpc)
            obj.count_nln(mpc);
            nk = count@mpe_branch_acc(obj, mpc);
        end
    end     %% methods
end         %% classdef
