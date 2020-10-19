classdef mpe_shunt_acp_nln < mpe_shunt_acp & mpe_wrapper_ac_nln

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpe_class = @mpe_shunt_acp;
    end

    methods
        function obj = mpe_shunt_acp_nln()
            obj@mpe_shunt_acp();
            obj.mpe_wrapper_ac_nln_init();
        end

        function build_params(obj, nm, mpc)
            build_params@mpe_shunt_acp(obj, nm, mpc);
            obj.build_nln_params(nm, mpc);
        end

        function nk = count(obj, mpc)
            obj.count_nln(mpc);
            nk = count@mpe_shunt_acp(obj, mpc);
        end
    end     %% methods
end         %% classdef
