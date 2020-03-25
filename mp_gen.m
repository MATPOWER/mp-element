classdef mp_gen < mp_element

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
        function obj = mp_gen(varargin)
            obj@mp_element(varargin{:});
            obj.name = 'gen';
            obj.mpc_field = 'gen';
            obj.np = 1;             %% this is a 1 port element
            obj.nz = 1;
        end

        function obj = add_states(obj, asm, mpc)
%             %% define constants
%             [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS] = idx_gen;

            ng = obj.nk;            %% number of gens
            asm.add_state(obj.name, ng);
        end

        function obj = build_params(obj, asm, mpc)
            %% define constants
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS] = idx_gen;

            %% incidence matrices
            nn = asm.getN('node');
            ng = obj.nk;
            IDs = mpc.gen(:, GEN_BUS);              %% bus IDs
            nidx = asm.node.data.ID2idx.bus(IDs);   %% node indexes
            obj.C = { sparse(nidx, 1:ng, 1, nn, ng) };

            nz = asm.getN('state');
            sidx = asm.state.data.ID2idx.(obj.name);    %% state indexes
            obj.D = { sparse(sidx, 1:ng, 1, nz, ng) };
        end

        %%-----  OPF methods  -----
        function build_gen_cost_params(obj, mpc, mpopt, pcost)
            %% find/prepare polynomial generator costs for active power
            [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
            ng = obj.nk;
            if nargin < 4
                [pcost qcost] = pqcost(mpc.gencost, ng);
            end
            cpg = [];
            cqg = [];
            ip0 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 1);   %% constant
            ip1 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 2);   %% linear
            ip2 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 3);   %% quadratic
            ip3 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) > 3);    %% cubic or greater
            if ~isempty(ip2) || ~isempty(ip1) || ~isempty(ip0)
                kpg = zeros(ng, 1);
                cpg = zeros(ng, 1);
                if ~isempty(ip2)
                    Qpg = zeros(ng, 1);
                    Qpg(ip2) = 2 * pcost(ip2, COST) * mpc.baseMVA^2;
                    cpg(ip2) = cpg(ip2) + pcost(ip2, COST+1) * mpc.baseMVA;
                    kpg(ip2) = kpg(ip2) + pcost(ip2, COST+2);
                else
                    Qpg = [];   %% no quadratic terms
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
                    'ip0', ip0, ...
                    'ip1', ip1, ...
                    'ip2', ip2, ...
                    'ip3', ip3, ...
                    'kpg', kpg, ...
                    'cpg', cpg, ...
                    'Qpg', Qpg ...
                );
            if ~isempty(ip3)
                error('mp_gen/build_gen_cost_params: polynomial generator costs greater than quadratic order not yet implemented');
            end

            %% find/prepare piecewise linear generator costs
            ipwl = find(mpc.gencost(:, MODEL) == PW_LINEAR);  %% piece-wise linear costs
            ny = size(ipwl, 1);   %% number of piece-wise linear cost vars
            nq = ng;    %% number of Qg variables
            [Ay, by] = makeAy(mpc.baseMVA, ng, mpc.gencost, 1, 1+ng, 1+ng+nq);
            obj.cost_pwl = struct('ny', ny, 'ipwl', ipwl, 'Ay', Ay, 'by', by);
        end
        
        function add_opf_vars(obj, asm, om, mpc, mpopt)
            %% collect/construct all generator cost parameters
            obj.build_gen_cost_params(mpc, mpopt);

            %% piecewise linear costs
            if obj.cost_pwl.ny
                om.add_var('y', obj.cost_pwl.ny);
            end
        end

        function add_opf_constraints(obj, asm, om, mpc, mpopt)
            %% piecewise linear costs
            if obj.cost_pwl.ny
                om.add_lin_constraint('ycon', obj.cost_pwl.Ay, [], obj.cost_pwl.by, {'Pg', 'Qg', 'y'});
            end
        end

        function add_opf_costs(obj, asm, om, mpc, mpopt)
            %% (quadratic) polynomial costs
            om.add_quad_cost('polPg', obj.cost_poly_p.Qpg, obj.cost_poly_p.cpg, obj.cost_poly_p.kpg, {'Pg'});

            %% piecewise linear costs
            if obj.cost_pwl.ny
                om.add_quad_cost('pwl', [], ones(obj.cost_pwl.ny, 1), 0, {'y'});
            end
        end
    end     %% methods
end         %% classdef
