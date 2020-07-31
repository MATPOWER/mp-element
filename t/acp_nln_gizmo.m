classdef acp_nln_gizmo < acp_gizmo & mpe_wrapper_ac_nln

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpe_class = @acp_gizmo;
    end

    methods
        function obj = acp_nln_gizmo()
            obj@acp_gizmo();
            obj.mpe_wrapper_ac_nln_init();
        end

        function build_params(obj, nm, mpc)
            build_params@acp_gizmo(obj, nm, mpc);
            obj.build_nln_params(nm, mpc);
        end

        function nk = count(obj, mpc)
            obj.count_nln(mpc);
            nk = count@acp_gizmo(obj, mpc);
        end
    end     %% methods
end         %% classdef
