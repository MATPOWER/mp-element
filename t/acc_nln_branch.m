classdef acc_nln_branch < acc_branch & ac_nln_wrapper

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpe_class = @acc_branch;
    end

    methods
        function obj = acc_nln_branch()
            obj@acc_branch();
            obj.ac_nln_wrapper_init();
        end

        function build_params(obj, asm, mpc)
            build_params@acc_branch(obj, asm, mpc);
            obj.build_nln_params(asm, mpc);
        end

        function nk = count(obj, mpc)
            obj.count_nln(mpc);
            nk = count@acc_branch(obj, mpc);
        end
    end     %% methods
end         %% classdef