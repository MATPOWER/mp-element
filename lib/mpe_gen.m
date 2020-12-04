classdef mpe_gen < mp_element

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        cost_pwl = struct();
        cost_poly_p = struct();
    end
    
    methods
        %% constructor
        function obj = mpe_gen()
            obj@mp_element();
            obj.name = 'gen';
            obj.dm_table = 'gen';
            obj.np = 1;             %% this is a 1 port element
            obj.nz = 1;
        end

        function obj = add_states(obj, nm, dm)
%             %% define constants
%             [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS] = idx_gen;

            ng = obj.nk;            %% number of gens
            nm.add_state(obj.name, ng);
        end

        function obj = build_params(obj, nm, dm)
            %% define constants
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS] = idx_gen;

            %% incidence matrices
            IDs = dm.mpc.gen(:, GEN_BUS);               %% bus IDs
            nidx = nm.node.data.ID2idx.bus(IDs);        %% node indexes
            sidx = nm.state.data.ID2idx.(obj.name);     %% state indexes
            obj.C = obj.incidence_matrix(nm.getN('node'), nidx);
            obj.D = obj.incidence_matrix(nm.getN('state'), sidx);
        end

        %%-----  OPF methods  -----
        function build_gen_cost_params(obj, mpc, mpopt, pcost)
            %% find/prepare polynomial generator costs for active power
            [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
            ng = obj.nk;
            if nargin < 4
                [pcost qcost] = pqcost(mpc.gencost, ng);
            end
            have_quad_cost = 0;
            kpg = []; cpg = []; Qpg = [];
            ip0 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 1);   %% constant
            ip1 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 2);   %% linear
            ip2 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 3);   %% quadratic
            ip3 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) > 3);    %% cubic or greater
            if ~isempty(ip2) || ~isempty(ip1) || ~isempty(ip0)
                have_quad_cost = 1;
                kpg = zeros(ng, 1);
                cpg = zeros(ng, 1);
                if ~isempty(ip2)
                    Qpg = zeros(ng, 1);
                    Qpg(ip2) = 2 * pcost(ip2, COST) * mpc.baseMVA^2;
                    cpg(ip2) = cpg(ip2) + pcost(ip2, COST+1) * mpc.baseMVA;
                    kpg(ip2) = kpg(ip2) + pcost(ip2, COST+2);
                end
                if ~isempty(ip1)
                    cpg(ip1) = cpg(ip1) + pcost(ip1, COST) * mpc.baseMVA;
                    kpg(ip1) = kpg(ip1) + pcost(ip1, COST+1);
                end
                if ~isempty(ip0)
                    kpg(ip0) = kpg(ip0) + pcost(ip0, COST);
                end
            end
            obj.cost_poly_p = struct( ...
                    'have_quad_cost', have_quad_cost, ...
                    'ip0', ip0, ...
                    'ip1', ip1, ...
                    'ip2', ip2, ...
                    'ip3', ip3, ...
                    'kpg', kpg, ...
                    'cpg', cpg, ...
                    'Qpg', Qpg ...
                );

            %% find/prepare piecewise linear generator costs
            ipwl = find(mpc.gencost(:, MODEL) == PW_LINEAR);  %% piece-wise linear costs
            ny = size(ipwl, 1);   %% number of piece-wise linear cost vars
            if isa(obj, 'mp_model_dc')
                nq = 0;    %% number of Qg variables
                q1 = [];
            else
                nq = ng;    %% number of Qg variables
                q1 = 1+ng;
            end
            [Ay, by] = makeAy(mpc.baseMVA, ng, mpc.gencost, 1, q1, 1+ng+nq);
            obj.cost_pwl = struct('ny', ny, 'ipwl', ipwl, 'Ay', Ay, 'by', by);
        end
        
        function add_opf_vars(obj, nm, om, dm, mpopt)
            %% collect/construct all generator cost parameters
            obj.build_gen_cost_params(dm.mpc, mpopt);

            %% piecewise linear costs
            if obj.cost_pwl.ny
                om.add_var('y', obj.cost_pwl.ny);
            end
        end

        function add_opf_costs(obj, nm, om, dm, mpopt)
            %% (quadratic) polynomial costs on Pg
            if obj.cost_poly_p.have_quad_cost
                om.add_quad_cost('polPg', obj.cost_poly_p.Qpg, obj.cost_poly_p.cpg, obj.cost_poly_p.kpg, {'Pg'});
            end

            %% (order 3 and higher) polynomial costs on Pg
            if ~isempty(obj.cost_poly_p.ip3)
                [pcost qcost] = pqcost(dm.mpc.gencost, obj.nk);
                cost_Pg = @(xx)opf_gen_cost_fcn(xx, dm.mpc.baseMVA, pcost, obj.cost_poly_p.ip3, mpopt);
                om.add_nln_cost('polPg', 1, cost_Pg, {'Pg'});
            end

            %% piecewise linear costs
            if obj.cost_pwl.ny
                om.add_quad_cost('pwl', [], ones(obj.cost_pwl.ny, 1), 0, {'y'});
            end
        end
    end     %% methods
end         %% classdef
