classdef acp_nln_gen < acp_gen & ac_nln_wrapper

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpe_class = @acp_gen;
    end
    
    methods
        function obj = acp_nln_gen()
            obj@acp_gen();
            obj.ac_nln_wrapper_init();
        end

        function build_params(obj, asm, mpc)
            build_params@acp_gen(obj, asm, mpc);
            obj.build_nln_params(asm, mpc);
        end

        function nk = count(obj, mpc)
            obj.count_nln(mpc);
            nk = count@acp_gen(obj, mpc);
        end
    end     %% methods
end         %% classdef
