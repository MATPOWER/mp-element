classdef mme_gen_opf < mme_gen

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'branch';
%     end
    
    methods
        function obj = add_vars(obj, mm, nm, dm, mpopt)
            nme = obj.network_model_element(nm);

            %% collect/construct all generator cost parameters
            nme.opf_build_gen_cost_params(dm);

            %% piecewise linear costs
            if nme.cost.pwl.n
                mm.add_var('y', nme.cost.pwl.n);
            end
        end

        function obj = add_costs(obj, mm, nm, dm, mpopt)
            nme = obj.network_model_element(nm);

            %% (quadratic) polynomial costs on Pg
            if nme.cost.poly_p.have_quad_cost
                mm.add_quad_cost('polPg', nme.cost.poly_p.Q, nme.cost.poly_p.c, nme.cost.poly_p.k, {'Pg'});
            end

            %% (order 3 and higher) polynomial costs on Pg
            if ~isempty(nme.cost.poly_p.i3)
                dme = obj.data_model_element(dm);
                cost_Pg = @(xx)poly_cost_fcn(nme, xx, dm.base_mva, dme.cost_pg.poly_coef, nme.cost.poly_p.i3);
                mm.add_nln_cost('polPg', 1, cost_Pg, {'Pg'});
            end

            %% piecewise linear costs
            if nme.cost.pwl.n
                mm.add_quad_cost('pwl', [], ones(nme.cost.pwl.n, 1), 0, {'y'});
            end
        end
    end     %% methods
end         %% classdef
