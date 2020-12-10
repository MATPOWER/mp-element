classdef nme_load_acc_nln < nme_load_acc & nme_wrapper_ac_nln

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpe_class = @nme_load_acc;
    end

    methods
        function obj = nme_load_acc_nln()
            obj@nme_load_acc();
            obj.nme_wrapper_ac_nln_init();
        end

        function build_params(obj, nm, dm)
            build_params@nme_load_acc(obj, nm, dm);
            obj.build_nln_params(nm, dm);
        end

        function nk = count(obj, dm)
            obj.mpe.count(dm);
            nk = count@nme_load_acc(obj, dm);
        end
    end     %% methods
end         %% classdef
