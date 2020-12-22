classdef nme_gen_ac < nme_gen% & mp_form_ac

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
        function obj = add_zvars(obj, nm, dm, idx)
            ng = obj.nk;
            dme = obj.data_model_element(dm);
            nm.add_var('zr', 'Pg', ng, dme.Pg0, dme.Pmin, dme.Pmax);
            nm.add_var('zi', 'Qg', ng, dme.Qg0, dme.Qmin, dme.Qmax);
        end

        function obj = build_params(obj, nm, dm)
            build_params@nme_gen(obj, nm, dm);      %% call parent
            ng = obj.nk;
            obj.N = -speye(ng);
        end

        %%-----  OPF methods  -----
        function opf_build_gen_cost_params(obj, dm)
            dme = obj.data_model_element(dm);
            obj.cost = dme.opf_build_gen_cost_params(dm, 0);
        end

        function obj = opf_add_constraints(obj, mm, nm, dm, mpopt)
            %% generator PQ capability curve constraints
            [Apqh, ubpqh, Apql, ubpql, Apqdata] = dm.gen_pq_capability_constraint();
            mm.add_lin_constraint('PQh', Apqh, [], ubpqh, {'Pg', 'Qg'});      %% npqh
            mm.add_lin_constraint('PQl', Apql, [], ubpql, {'Pg', 'Qg'});      %% npql
            mm.userdata.Apqdata = Apqdata;

            %% dispatchable load constant power factor constraint
            [Avl, lvl, uvl] = dm.elm_by_name('gen').disp_load_constant_pf_constraint(dm);
            if ~isempty(Avl)
                mm.add_lin_constraint('vl',  Avl, lvl, uvl,   {'Pg', 'Qg'});    %% nvl
            end

            %% piecewise linear costs
            if obj.cost.pwl.n
                mm.add_lin_constraint('ycon', obj.cost.pwl.A, [], obj.cost.pwl.b, {'Pg', 'Qg', 'y'});
            end

            %% call parent
            opf_add_constraints@nme_gen(obj, mm, nm, dm, mpopt);
        end

        function obj = opf_add_costs(obj, mm, nm, dm, mpopt)
            %% call parent
            opf_add_costs@nme_gen(obj, mm, nm, dm, mpopt);

            %% costs on reactive dispatch
            if ~isempty(obj.cost.poly_q)
                %% (quadratic) polynomial costs on Qg
                if obj.cost.poly_q.have_quad_cost
                    mm.add_quad_cost('polQg', obj.cost.poly_q.Q, obj.cost.poly_q.c, obj.cost.poly_q.k, {'Qg'});
                end

                %% (order 3 and higher) polynomial costs on Qg
                if ~isempty(obj.cost.poly_q.i3)
                    dme = obj.data_model_element(dm);
                    cost_Qg = @(xx)opf_gen_cost_fcn(xx, obj.cost.baseMVA, dme.qcost, obj.cost.poly_q.i3);
                    mm.add_nln_cost('polQg', 1, cost_Qg, {'Qg'});
                end
            end
        end
    end     %% methods
end         %% classdef
