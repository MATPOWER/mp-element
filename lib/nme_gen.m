classdef nme_gen < nm_element

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        cost        %% struct for cost parameters with fields
                    %%  .poly_p - polynomial costs for active power,
                    %%            struct with fields:
                    %%      .have_quad_cost
                    %%      .i0, .i1, .i2, .i3
                    %%      .k, .c, .Q
                    %%  .poly_q - polynomial costs for reactive power,
                    %%            same struct as .poly_p
                    %%  .pwl - piecewise linear costs for actve & reactive
                    %%         struct with fields:
                    %%      .n, .i, .A, .b
    end
    
    methods
        %% constructor
        function obj = nme_gen()
            obj@nm_element();
            obj.name = 'gen';
            obj.np = 1;             %% this is a 1 port element
            obj.nz = 1;
        end

        function obj = add_states(obj, nm, dm)
            ng = obj.nk;            %% number of gens
            nm.add_state(obj.name, ng);
        end

        function idx = node_indices(obj, nm, dm, dme)
            bus_dme = dm.elements.bus;
            nidx = nm.get_node_idx('bus');  %% node indices for 'bus'
            b = dme.bus(dme.on);
            bidx = bus_dme.i2on(b);         %% online bus indices for gens
            idx = nidx(bidx);               %% node indices for gens
        end

        function obj = build_params(obj, nm, dm)
            dme = obj.data_model_element(dm);

            %% incidence matrices
            nidx = obj.node_indices(nm, dm, dme);
            sidx = nm.get_state_idx(obj.name);  %% state indices for 'gen'
            obj.C = obj.incidence_matrix(nm.getN('node'), nidx);
            obj.D = obj.incidence_matrix(nm.getN('state'), sidx);
        end

        %%-----  OPF methods  -----
        function obj = opf_add_vars(obj, mm, nm, dm, mpopt)
            %% collect/construct all generator cost parameters
            obj.opf_build_gen_cost_params(dm);

            %% piecewise linear costs
            if obj.cost.pwl.n
                mm.add_var('y', obj.cost.pwl.n);
            end
        end

        function obj = opf_add_costs(obj, mm, nm, dm, mpopt)
            %% (quadratic) polynomial costs on Pg
            if obj.cost.poly_p.have_quad_cost
                mm.add_quad_cost('polPg', obj.cost.poly_p.Q, obj.cost.poly_p.c, obj.cost.poly_p.k, {'Pg'});
            end

            %% (order 3 and higher) polynomial costs on Pg
            if ~isempty(obj.cost.poly_p.i3)
                dme = obj.data_model_element(dm);
                cost_Pg = @(xx)poly_cost_fcn(obj, xx, dm.base_mva, dme.cost_pg.poly_coef, obj.cost.poly_p.i3);
                mm.add_nln_cost('polPg', 1, cost_Pg, {'Pg'});
            end

            %% piecewise linear costs
            if obj.cost.pwl.n
                mm.add_quad_cost('pwl', [], ones(obj.cost.pwl.n, 1), 0, {'y'});
            end
        end

        function [f, df, d2f] = poly_cost_fcn(obj, xx, x_scale, cost, idx)
            x = xx{1}(idx) * x_scale;
            n = length(xx{1});

            %%----- evaluate cost function -----
            f = sum( eval_poly_fcn(cost(idx, :), x) );

            %%----- evaluate cost gradient -----
            if nargout > 1
                %% coefficients of 1st derivative
                cp = diff_poly_fcn(cost(idx, :));
                df = zeros(n, 1);
                df(idx) = x_scale * eval_poly_fcn(cp, x);

                %% ---- evaluate cost Hessian -----
                if nargout > 2
                    %% coefficients of 2nd derivative
                    cpp = diff_poly_fcn(cp);
                    d2f = sparse(idx, idx, x_scale^2 * eval_poly_fcn(cpp, x), n, n);
                end
            end
        end
    end     %% methods
end         %% classdef

function c = diff_poly_fcn(c)
    n = size(c, 2);     %% number of coefficients (cols in c)
    if n >= 2
        c = c(:, 2:n);
    else
        c = zeros(size(c, 1), 1);
    end
    for k = 2:n-1
        c(:, k) = k * c(:, k);
    end
end

function f = eval_poly_fcn(c, x)
    if isempty(c)
        f = zeros(size(x));
    else
        f = c(:, 1);        %% constant term
        for k = 2:size(c, 2)
            f = f + c(:, k) .* x .^ (k-1);
        end
    end
end
