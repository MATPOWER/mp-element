classdef dme_gen_opf < dme_gen
%DME_GEN_OPF  MATPOWER data model class for gen data

%   MATPOWER
%   Copyright (c) 1996-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        cost_pg     %% cost parameter table for active power generation, all gens
        cost_qg     %% cost parameter table for reactive power generation, all gens
        pwl1        %% indices of single-block piecewise linear costs, all gens
                    %% (automatically converted to linear cost)
    end     %% properties

    methods
        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dme_gen(obj), ...
                {'mu_pg_lb', 'mu_pg_ub', 'mu_qg_lb', 'mu_qg_ub'});
        end

        function vars = export_vars(obj, task)
            vars = horzcat( export_vars@dme_gen(obj), ...
                {'vm_setpoint', 'mu_pg_lb', 'mu_pg_ub', 'mu_qg_lb', 'mu_qg_ub'} );
        end

        function TorF = have_cost(obj)
            TorF = 1;
        end

        function cost = build_cost_params(obj, dm, dc)
            base_mva = dm.base_mva;

            poly_p = obj.gen_cost_poly_params(base_mva, obj.cost_pg(obj.on, :));
            if dc || isempty(obj.cost_qg)
                pwl = obj.gen_cost_pwl_params(base_mva, obj.cost_pg(obj.on, :), obj.n, dc);
                poly_q = [];
            else
                poly_q = obj.gen_cost_poly_params(base_mva, obj.cost_qg(obj.on, :));

                %% expand cost params as needed
                polyNp = max(obj.cost_pg.poly_n);
                polyNq = max(obj.cost_qg.poly_n);
                pwlNp  = max(obj.cost_pg.pwl_n);
                pwlNq  = max(obj.cost_qg.pwl_n);
                polyN = max(polyNp, polyNq);
                pwlN = max(pwlNp, pwlNq);
                if size(obj.cost_pg.poly_coef, 2) < polyN
                    obj.cost_pg.poly_coef(end, polyN) = 0;
                end
                if size(obj.cost_qg.poly_coef, 2) < polyN
                    obj.cost_qg.poly_coef(end, polyN) = 0;
                end
                if size(obj.cost_pg.pwl_cost, 2) < pwlN
                    obj.cost_pg.pwl_qty(end, pwlN) = 0;
                    obj.cost_pg.pwl_cost(end, pwlN) = 0;
                end
                if size(obj.cost_qg.pwl_cost, 2) < pwlN
                    obj.cost_qg.pwl_qty(end, pwlN) = 0;
                    obj.cost_qg.pwl_cost(end, pwlN) = 0;
                end

                %% stack pg & qg cost params first
                cost = [obj.cost_pg; obj.cost_qg];
                pwl = obj.gen_cost_pwl_params(base_mva, cost, obj.n, dc);
            end

            cost = struct( ...
                    'poly_p',   poly_p, ...
                    'poly_q',   poly_q, ...
                    'pwl',      pwl ...
                );
        end

        function p = gen_cost_pwl_params(obj, base_mva, cost, ng, dc)
            ipwl = find(cost.pwl_n);    %% piece-wise linear costs
            ny = size(ipwl, 1);     %% number of piece-wise linear cost vars
            if dc
                nq = 0;     %% number of qg variables
                q1 = [];
            else
                nq = ng;    %% number of qg variables
                q1 = 1+ng;
            end

            %% from makeAy()
            ybas = 1+ng+nq;
            if ny == 0
                Ay = sparse([], [], [], 0, ybas+ny-1, 0);
                by = [];
            else
                %% if p(i),p(i+1),c(i),c(i+1) define one of the cost segments,
                %% then the corresponding constraint on pg (or qg) and Y is
                %%                                             c(i+1) - c(i)
                %%  Y   >=   c(i) + m * (pg - p(i)),      m = ---------------
                %%                                             p(i+1) - p(i)
                %%
                %% this becomes   m * pg - Y   <=   m*p(i) - c(i)

                %% form constraint matrix
                m = sum(cost.pwl_n(ipwl));  %% total number of cost points
                Ay = sparse([], [], [], m-ny, ybas+ny-1, 2*(m-ny)); 
                by = [];
                k = 1;
                j = 1;
                for i=ipwl'
                    ns = cost.pwl_n(i); %% # of cost points; segments = ns-1
                    p = cost.pwl_qty(i, 1:ns) / base_mva;
                    c = cost.pwl_cost(i, 1:ns);
                    m = diff(c) ./ diff(p);         %% slopes for pg (or qg)
                    if any(diff(p) == 0)
                        fprintf('dme_gen/gen_cost_pwl_params: bad qty data in row %i of cost matrix\n', i);
                    end
                    b = m .* p(1:ns-1) - c(1:ns-1); %% and rhs
                    by = [by;  b'];
                    if i > ng
                        sidx = q1 + (i-ng) - 1;     %% for qg cost
                    else
                        sidx = i;                   %% for pg cost
                    end
                    Ay(k:k+ns-2, sidx) = m';
                    Ay(k:k+ns-2, ybas+j-1) = -ones(ns-1,1);
                    k = k + ns - 1;
                    j = j + 1;
                end
            end

            p = struct('n', ny, 'i', ipwl, 'A', Ay, 'b', by);
        end

        function p = gen_cost_poly_params(obj, base_mva, cost)
            ng = size(cost, 1);
            have_quad_cost = 0;
            kg = []; cg = []; Qg = [];
            i0 = find(cost.poly_n == 1);    %% constant
            i1 = find(cost.poly_n == 2);    %% linear
            i2 = find(cost.poly_n == 3);    %% quadratic
            i3 = find(cost.poly_n > 3);     %% cubic or greater
            if ~isempty(i2) || ~isempty(i1) || ~isempty(i0)
                have_quad_cost = 1;
                kg = zeros(ng, 1);
                cg = zeros(ng, 1);
                if ~isempty(i2)
                    Qg = zeros(ng, 1);
                    Qg(i2) = 2 * cost.poly_coef(i2, 3) * base_mva^2;
                    cg(i2) = cg(i2) + cost.poly_coef(i2, 2) * base_mva;
                    kg(i2) = kg(i2) + cost.poly_coef(i2, 1);
                end
                if ~isempty(i1)
                    cg(i1) = cg(i1) + cost.poly_coef(i1, 2) * base_mva;
                    kg(i1) = kg(i1) + cost.poly_coef(i1, 1);
                end
                if ~isempty(i0)
                    kg(i0) = kg(i0) + cost.poly_coef(i0, 1);
                end
            end
            p = struct( ...
                    'have_quad_cost', have_quad_cost, ...
                    'i0', i0, ...
                    'i1', i1, ...
                    'i2', i2, ...
                    'i3', i3, ...
                    'k', kg, ...
                    'c', cg, ...
                    'Q', Qg ...
                );
        end

        function maxgc = max_pwl_gencost(obj)
            maxgc = max(max(obj.cost_pg.pwl_cost));
            if ~isempty(obj.cost_qg) && ~isempty(obj.cost_qg.pwl_cost)
                maxgc = max(maxgc, max(max(obj.cost_qg.pwl_cost)));
            end
        end
    end     %% methods
end         %% classdef