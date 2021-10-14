classdef mme_gen_opf_ac < mme_gen_opf

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

            %% generator PQ capability curve constraints
            [Apqh, ubpqh, Apql, ubpql, Apqdata] = ...
                dm.elements.gen.pq_capability_constraint(dm.base_mva);
            mm.add_lin_constraint('PQh', Apqh, [], ubpqh, {'Pg', 'Qg'});      %% npqh
            mm.add_lin_constraint('PQl', Apql, [], ubpql, {'Pg', 'Qg'});      %% npql
            mm.userdata.Apqdata = Apqdata;

            %% dispatchable load constant power factor constraint
            [Avl, lvl, uvl] = dm.elements.gen.disp_load_constant_pf_constraint(dm);
            if ~isempty(Avl)
                mm.add_lin_constraint('vl',  Avl, lvl, uvl,   {'Pg', 'Qg'});    %% nvl
            end

            %% piecewise linear costs
            if nme.cost.pwl.n
                mm.add_lin_constraint('ycon', nme.cost.pwl.A, [], nme.cost.pwl.b, {'Pg', 'Qg', 'y'});
            end

            %% call parent
            add_constraints@mme_gen_opf(obj, mm, nm, dm, mpopt);
        end

        function obj = add_costs(obj, mm, nm, dm, mpopt)
            nme = obj.network_model_element(nm);

            %% call parent
            add_costs@mme_gen_opf(obj, mm, nm, dm, mpopt);

            %% costs on reactive dispatch
            if ~isempty(nme.cost.poly_q)
                %% (quadratic) polynomial costs on Qg
                if nme.cost.poly_q.have_quad_cost
                    mm.add_quad_cost('polQg', nme.cost.poly_q.Q, nme.cost.poly_q.c, nme.cost.poly_q.k, {'Qg'});
                end

                %% (order 3 and higher) polynomial costs on Qg
                if ~isempty(nme.cost.poly_q.i3)
                    dme = obj.data_model_element(dm);
                    cost_Qg = @(xx)poly_cost_fcn(nme, xx, dm.base_mva, dme.cost_qg.poly_coef, nme.cost.poly_q.i3);
                    mm.add_nln_cost('polQg', 1, cost_Qg, {'Qg'});
                end
            end
        end
    end     %% methods
end         %% classdef
