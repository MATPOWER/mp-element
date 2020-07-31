classdef acp_nln_branch < acp_branch & ac_nln_wrapper

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpe_class = @acp_branch;
    end

    methods
        function obj = acp_nln_branch()
            obj@acp_branch();
            obj.ac_nln_wrapper_init();
        end

        function build_params(obj, nm, mpc)
            build_params@acp_branch(obj, nm, mpc);
            obj.build_nln_params(nm, mpc);
        end

        function nk = count(obj, mpc)
            obj.count_nln(mpc);
            nk = count@acp_branch(obj, mpc);
        end
    end     %% methods
end         %% classdef
