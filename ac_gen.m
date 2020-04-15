classdef ac_gen < mp_gen% & ac_model

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        cost_poly_q = struct();
    end
    
    methods
        function obj = add_zvars(obj, asm, mpc, idx)
            %% define constants
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            ng = obj.nk;
            Pg   = mpc.gen(:, PG) / mpc.baseMVA;
            Pmin = mpc.gen(:, PMIN) / mpc.baseMVA;
            Pmax = mpc.gen(:, PMAX) / mpc.baseMVA;
            Qg   = mpc.gen(:, QG) / mpc.baseMVA;
            Qmin = mpc.gen(:, QMIN) / mpc.baseMVA;
            Qmax = mpc.gen(:, QMAX) / mpc.baseMVA;
            asm.add_var('zr', 'Pg', ng, Pg, Pmin, Pmax);
            asm.add_var('zi', 'Qg', ng, Qg, Qmin, Qmax);
        end

        function obj = build_params(obj, asm, mpc)
            build_params@mp_gen(obj, asm, mpc);     %% call parent
            ng = obj.nk;
            obj.N = -speye(ng);
        end

        function build_gen_cost_params(obj, mpc, mpopt, pcost, qcost)
            [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
            ng = obj.nk;
            if nargin < 5
                [pcost qcost] = pqcost(mpc.gencost, ng);
            end

            %% call parent (with pcost) to handle active power costs
            build_gen_cost_params@mp_gen(obj, mpc, mpopt, pcost);
            
            %% find/prepare polynomial generator costs for reactive power
            obj.cost_poly_q.have_quad_cost = 0;
            obj.cost_poly_q.iq3 = [];
            if ~isempty(qcost)
                have_quad_cost = 0;
                kqg = []; cqg = []; Qqg = [];
                iq0 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) == 1);   %% constant
                iq1 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) == 2);   %% linear
                iq2 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) == 3);   %% quadratic
                iq3 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) > 3);    %% cubic or greater
                if ~isempty(iq2) || ~isempty(iq1) || ~isempty(iq0)
                    have_quad_cost = 1;
                    kqg = zeros(ng, 1);
                    cqg = zeros(ng, 1);
                    if ~isempty(iq2)
                        Qqg = zeros(ng, 1);
                        Qqg(iq2) = 2 * qcost(iq2, COST) * mpc.baseMVA^2;
                        cqg(iq2) = cqg(iq2) + qcost(iq2, COST+1) * mpc.baseMVA;
                        kqg(iq2) = kqg(iq2) + qcost(iq2, COST+2);
                    end
                    if ~isempty(iq1)
                        cqg(iq1) = cqg(iq1) + qcost(iq1, COST) * mpc.baseMVA;
                        kqg(iq1) = kqg(iq1) + qcost(iq1, COST+1);
                    end
                    if ~isempty(iq0)
                        kqg(iq0) = kqg(iq0) + qcost(iq0, COST);
                    end
                end
                obj.cost_poly_q = struct( ...
                        'have_quad_cost', have_quad_cost, ...
                        'iq0', iq0, ...
                        'iq1', iq1, ...
                        'iq2', iq2, ...
                        'iq3', iq3, ...
                        'kqg', kqg, ...
                        'cqg', cqg, ...
                        'Qqg', Qqg ...
                    );
            end
        end

        function add_opf_constraints(obj, asm, om, mpc, mpopt)
            %% generator PQ capability curve constraints
            [Apqh, ubpqh, Apql, ubpql, Apqdata] = makeApq(mpc.baseMVA, mpc.gen);
            om.add_lin_constraint('PQh', Apqh, [], ubpqh, {'Pg', 'Qg'});      %% npqh
            om.add_lin_constraint('PQl', Apql, [], ubpql, {'Pg', 'Qg'});      %% npql
            om.userdata.Apqdata = Apqdata;

            %% dispatchable load constant power factor constraint
            [Avl, lvl, uvl]  = makeAvl(mpc);
            if ~isempty(Avl)
                om.add_lin_constraint('vl',  Avl, lvl, uvl,   {'Pg', 'Qg'});    %% nvl
            end

            %% piecewise linear costs
            if obj.cost_pwl.ny
                om.add_lin_constraint('ycon', obj.cost_pwl.Ay, [], obj.cost_pwl.by, {'Pg', 'Qg', 'y'});
            end

            %% call parent
            add_opf_constraints@mp_gen(obj, asm, om, mpc, mpopt);
        end

        function add_opf_costs(obj, asm, om, mpc, mpopt)
            %% call parent
            add_opf_costs@mp_gen(obj, asm, om, mpc, mpopt);

            %% (quadratic) polynomial costs on Qg
            if obj.cost_poly_q.have_quad_cost
                om.add_quad_cost('polQg', obj.cost_poly_q.Qqg, obj.cost_poly_q.cqg, obj.cost_poly_q.kqg, {'Qg'});
            end

            %% (order 3 and higher) polynomial costs on Qg
            if ~isempty(obj.cost_poly_q.iq3)
                [pcost qcost] = pqcost(mpc.gencost, obj.nk);
                cost_Qg = @(xx)opf_gen_cost_fcn(xx, mpc.baseMVA, qcost, obj.cost_poly_q.iq3, mpopt);
                om.add_nln_cost('polQg', 1, cost_Qg, {'Qg'});
            end
        end
    end     %% methods
end         %% classdef
