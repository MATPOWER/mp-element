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

        function obj = build_params(obj, nm, dm)
            dme = obj.data_model_element(dm);

            %% incidence matrices
            nidx = nm.get_node_idx('bus');      %% node indices for 'bus'
            sidx = nm.get_state_idx(obj.name);  %% state indices for 'gen'
            idx = nidx(dme.bus);                %% node indices for gens
            obj.C = obj.incidence_matrix(nm.getN('node'), idx);
            obj.D = obj.incidence_matrix(nm.getN('state'), sidx);
        end

        %%-----  OPF methods  -----
        function add_opf_vars(obj, mm, nm, dm, mpopt)
            %% collect/construct all generator cost parameters
            obj.build_gen_cost_params(dm);

            %% piecewise linear costs
            if obj.cost.pwl.n
                mm.add_var('y', obj.cost.pwl.n);
            end
        end

        function add_opf_costs(obj, mm, nm, dm, mpopt)
            %% (quadratic) polynomial costs on Pg
            if obj.cost.poly_p.have_quad_cost
                mm.add_quad_cost('polPg', obj.cost.poly_p.Q, obj.cost.poly_p.c, obj.cost.poly_p.k, {'Pg'});
            end

            %% (order 3 and higher) polynomial costs on Pg
            if ~isempty(obj.cost.poly_p.i3)
                dme = obj.data_model_element(dm);
                cost_Pg = @(xx)opf_gen_cost_fcn(xx, obj.cost.baseMVA, dme.pcost, obj.cost.poly_p.i3);
                mm.add_nln_cost('polPg', 1, cost_Pg, {'Pg'});
            end

            %% piecewise linear costs
            if obj.cost.pwl.n
                mm.add_quad_cost('pwl', [], ones(obj.cost.pwl.n, 1), 0, {'y'});
            end
        end
    end     %% methods
end         %% classdef
