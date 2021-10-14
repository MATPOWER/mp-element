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
            nme = obj.network_model_element(nm);

            %% piecewise linear costs
            if nme.cost.pwl.n
                mm.add_lin_constraint('ycon', nme.cost.pwl.A, [], nme.cost.pwl.b, {'Pg', 'y'});
            end

            %% call parent
            add_constraints@mme_gen_opf(obj, mm, nm, dm, mpopt);
        end
    end     %% methods
end         %% classdef
