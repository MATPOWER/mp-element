classdef mme_reserve_gen < mp.mm_element

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
        function name = name(obj)
            name = 'reserve_gen';
        end

        function obj = add_vars(obj, mm, nm, dm, mpopt)
            dme = obj.data_model_element(dm);
            mm.add_var('R', dme.n, 0, 0, dme.r_ub);
        end

        function obj = add_costs(obj, mm, nm, dm, mpopt)
            dme = obj.data_model_element(dm);
            c = dme.tab.cost(dme.on) * dm.base_mva;    %% p.u. cost coeffs
            mm.add_quad_cost('Rcost', [], c, 0, {'R'});
        end

        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            dme = obj.data_model_element(dm);
            gen_dme = dm.elements.gen;
            ng = gen_dme.n; %% number of online generators
            ngr = dme.n;    %% number of online reserve gens
            igr = dme.gen;  %% index of online gen for each online reserve gen

            Ar = [sparse(1:ngr, igr, 1, ngr, ng) speye(ngr)];
            ur = gen_dme.pg_ub(igr);

            mm.add_lin_constraint('Pg_plus_R', Ar, [], ur, {'Pg', 'R'});
        end

        function obj = data_model_update(obj, mm, nm, dm, mpopt)
            dme = obj.data_model_element(dm);
            rz_dme = dm.elements.reserve_zone;

            %% get reserve allocations, bounds, shadow prices, cost
            [vv, ll] = mm.get_idx('var', 'lin');
            r = mm.soln.x(vv.i1.R:vv.iN.R) * dm.base_mva;
            r_lb = 0;
            r_ub = dme.r_ub * dm.base_mva;
            zprc = zeros(rz_dme.nr, 1);
            zprc(rz_dme.on) = mm.soln.lambda.mu_l(ll.i1.Rreq:ll.iN.Rreq) / dm.base_mva;
            prc = rz_dme.tab.zones(:, dme.tab.gen)' * zprc;
            mu_lb = mm.soln.lambda.lower(vv.i1.R:vv.iN.R) / dm.base_mva;
            mu_ub = mm.soln.lambda.upper(vv.i1.R:vv.iN.R) / dm.base_mva;
            mu_pg_ub = mm.soln.lambda.mu_u(ll.i1.Pg_plus_R:ll.iN.Pg_plus_R) / dm.base_mva;
            total_cost = mm.eval_quad_cost(mm.soln.x, 'Rcost');

            %% update in the data model
            dme.tab.r(dme.on) = r;
            dme.tab.r_lb(dme.on) = r_lb;
            dme.tab.r_ub(dme.on) = r_ub;
            dme.tab.prc = prc;
            dme.tab.mu_lb(dme.on) = mu_lb;
            dme.tab.mu_ub(dme.on) = mu_ub;
            dme.tab.mu_pg_ub(dme.on) = mu_pg_ub;
            dme.tab.total_cost(dme.on) = total_cost;
        end
    end     %% methods
end         %% classdef
