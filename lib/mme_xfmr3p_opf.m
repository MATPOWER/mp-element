classdef mme_xfmr3p_opf < mme_xfmr3p

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end
    
    methods
        function x0 = opf_interior_x0(obj, mm, nm, dm, x0)
        end

        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
            nme = obj.network_model_element(nm);
            nme.pf_data_model_update(mm, nm, dm, mpopt);
        end
    end     %% methods
end         %% classdef
