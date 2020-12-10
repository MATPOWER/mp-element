classdef nme_shunt_dc < nme_shunt & mp_form_dc

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
            build_params@nme_shunt(obj, nm, dm);   %% call parent

            dme = obj.data_model_element(dm);
            obj.p = dme.Gs(dme.on);     %% vector of shunt conductances
        end
    end     %% methods
end         %% classdef
