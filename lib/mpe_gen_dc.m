classdef mpe_gen_dc < mpe_gen & mp_model_dc

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'gen';
%     end
    
    methods
        function obj = add_zvars(obj, nm, dm, idx)
            ng = obj.nk;
            dme = obj.data_model_element(dm);
            nm.add_var('z', 'Pg', ng, dme.Pg0, dme.Pmin, dme.Pmax);
        end

        function obj = build_params(obj, nm, dm)
            build_params@mpe_gen(obj, nm, dm);      %% call parent
            ng = obj.nk;
            obj.K = -speye(ng);
        end

        function add_opf_constraints(obj, nm, om, dm, mpopt)
            %% piecewise linear costs
            if obj.cost.pwl.n
                om.add_lin_constraint('ycon', obj.cost.pwl.A, [], obj.cost.pwl.b, {'Pg', 'y'});
            end

            %% call parent
            add_opf_constraints@mpe_gen(obj, nm, om, dm, mpopt);
        end
    end     %% methods
end         %% classdef
