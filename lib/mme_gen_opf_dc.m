classdef mme_gen_opf_dc < mme_gen_opf

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'gen';
%     end
    
    methods
        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %% piecewise linear costs
            if obj.cost.pwl.n
                mm.add_lin_constraint('ycon', obj.cost.pwl.A, [], obj.cost.pwl.b, {'Pg', 'y'});
            end

            %% call parent
            add_constraints@mme_gen_opf(obj, mm, nm, dm, mpopt);
        end

        function opf_build_gen_cost_params(obj, dm)
            dme = obj.data_model_element(dm);
            obj.cost = dme.opf_build_gen_cost_params(dm, 1);
        end
    end     %% methods
end         %% classdef
