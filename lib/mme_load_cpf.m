classdef mme_load_cpf < mm_element

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'load';
%     end
    
    methods
        %% constructor
        function obj = mme_load_cpf()
            obj@mm_element();
            obj.name = 'load';
        end

        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            ad = mm.aux_data;
            dme = obj.data_model_element(dm);
            dm = dme.parameterized(dm, ad.dmb, ad.dmt, mm.soln.x(end));
        end
    end     %% methods
end         %% classdef
