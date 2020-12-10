classdef nme_load_dc < nme_load & mp_model_dc

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = build_params(obj, nm, dm)
            build_params@nme_load(obj, nm, dm);     %% call parent

            dme = obj.data_model_element(dm);
            obj.p = dme.Pd(dme.on);     %% vector of active power demand
        end
    end     %% methods
end         %% classdef
