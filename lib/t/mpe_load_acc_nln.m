classdef mpe_load_acc_nln < mpe_load_acc & mpe_wrapper_ac_nln

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpe_class = @mpe_load_acc;
    end

    methods
        function obj = mpe_load_acc_nln()
            obj@mpe_load_acc();
            obj.mpe_wrapper_ac_nln_init();
        end

        function build_params(obj, nm, dm)
            build_params@mpe_load_acc(obj, nm, dm);
            obj.build_nln_params(nm, dm);
        end

        function nk = count(obj, dm)
            obj.count_nln(dm);
            nk = count@mpe_load_acc(obj, dm);
        end
    end     %% methods
end         %% classdef
