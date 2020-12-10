classdef nme_gen_acp_nln < nme_gen_acp & nme_wrapper_ac_nln

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpe_class = @nme_gen_acp;
    end
    
    methods
        function obj = nme_gen_acp_nln()
            obj@nme_gen_acp();
            obj.nme_wrapper_ac_nln_init();
        end

        function build_params(obj, nm, dm)
            build_params@nme_gen_acp(obj, nm, dm);
            obj.build_nln_params(nm, dm);
        end

        function nk = count(obj, dm)
            obj.mpe.count(dm);
            nk = count@nme_gen_acp(obj, dm);
        end
    end     %% methods
end         %% classdef
