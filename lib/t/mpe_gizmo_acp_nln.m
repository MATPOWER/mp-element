classdef mpe_gizmo_acp_nln < mpe_gizmo_acp & mpe_wrapper_ac_nln

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpe_class = @mpe_gizmo_acp;
    end

    methods
        function obj = mpe_gizmo_acp_nln()
            obj@mpe_gizmo_acp();
            obj.mpe_wrapper_ac_nln_init();
        end

        function build_params(obj, nm, dm)
            build_params@mpe_gizmo_acp(obj, nm, dm);
            obj.build_nln_params(nm, dm);
        end

        function nk = count(obj, dm)
            obj.count_nln(dm);
            nk = count@mpe_gizmo_acp(obj, dm);
        end
    end     %% methods
end         %% classdef
