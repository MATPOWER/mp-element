classdef mpe_bus_dc < mpe_bus & mp_model_dc

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'bus';
%     end
    
    methods
        function obj = add_vvars(obj, nm, dm, idx)
            dme = obj.data_model_element(dm);
            nb = obj.nk;

            %% prepare angle bounds for ref buses
            Vamin = -Inf(nb, 1);
            Vamax =  Inf(nb, 1);
            k = find(dme.isref(dme.on));
            Vamin(k) = dme.Va0(k);
            Vamax(k) = dme.Va0(k);

            nm.add_var('va', 'Va', nb, dme.Va0, Vamin, Vamax);
        end
    end     %% methods
end         %% classdef
