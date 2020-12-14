classdef dme_gen_mpc2 < dme_gen & dm_format_mpc2
%DME_GEN_MPC2  MATPOWER data model gen table for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        pwl1        %% indices of single-block piecewise linear costs
    end     %% properties

    methods
        %% constructor
        function obj = dme_gen_mpc2()
            obj@dme_gen();      %% call parent constructor

            %% define named indices into data matrices
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS] = idx_gen;
            obj.st_col = GEN_STATUS;
        end

        function obj = initialize(obj, dm)
            initialize@dme_gen(obj, dm);    %% call parent

            %% define named indices into data matrices
            [GEN_BUS] = idx_gen;

            %% get bus mapping info
            b2i = dm.elm_by_name('bus').ID2i;   %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            tab = obj.get_table(dm);
            obj.bus = b2i(tab(:, GEN_BUS));
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elm_by_name('bus').status;  %% bus status

            %% update status of gens at isolated/offline buses
            obj.status = obj.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@dme_gen(obj, dm);
        end

        function obj = build_params(obj, dm)
            %% define named indices into data matrices
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
            [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
            baseMVA = dm.mpc.baseMVA;

            gen = obj.get_table(dm);

            %% get generator parameters
            obj.Pg0  = gen(obj.on, PG) / baseMVA;
            obj.Pmin = gen(obj.on, PMIN) / baseMVA;
            obj.Pmax = gen(obj.on, PMAX) / baseMVA;
            obj.Qg0  = gen(obj.on, QG) / baseMVA;
            obj.Qmin = gen(obj.on, QMIN) / baseMVA;
            obj.Qmax = gen(obj.on, QMAX) / baseMVA;
            obj.Vg   = gen(obj.on, VG);

            %% get gen cost parameters
            if isfield(dm.mpc, 'gencost') && ~isempty(dm.mpc.gencost)
                gencost = dm.mpc.gencost;

                %% convert single-block piecewise-linear costs into linear polynomial cost
                pwl1 = find(gencost(:, MODEL) == PW_LINEAR & gencost(:, NCOST) == 2);
                % p1 = [];
                if ~isempty(pwl1)
                    x0 = gencost(pwl1, COST);
                    y0 = gencost(pwl1, COST+1);
                    x1 = gencost(pwl1, COST+2);
                    y1 = gencost(pwl1, COST+3);
                    m = (y1 - y0) ./ (x1 - x0);
                    b = y0 - m .* x0;
                    gencost(pwl1, MODEL) = POLYNOMIAL;
                    gencost(pwl1, NCOST) = 2;
                    gencost(pwl1, COST:COST+1) = [m b];
                    obj.pwl1 = pwl1;
                end

                [pcost, qcost] = pqcost(gencost, obj.nr);
                obj.pcost = pcost(obj.on, :);
                if isempty(qcost)
                    obj.qcost = qcost;
                else
                    obj.qcost = qcost(obj.on, :);
                end
            end
        end

        %%-----  OPF methods  -----
        function cost = build_gen_cost_params(obj, dm, dc)
            mpc = dm.mpc;
            baseMVA = dm.mpc.baseMVA;

            poly_p = obj.gen_cost_poly_params(baseMVA, obj.pcost);
            if dc || isempty(obj.qcost)
                poly_q = [];
                pwl = obj.gen_cost_pwl_params(baseMVA, obj.pcost, obj.n, dc);
            else
                poly_q = obj.gen_cost_poly_params(baseMVA, obj.qcost);
                pwl = obj.gen_cost_pwl_params(baseMVA, ...
                            [obj.pcost; obj.qcost], obj.n, dc);
            end

            cost = struct( ...
                    'baseMVA',  mpc.baseMVA, ...
                    'poly_p',   poly_p, ...
                    'poly_q',   poly_q, ...
                    'pwl',      pwl ...
                );
        end

        function p = gen_cost_pwl_params(obj, baseMVA, cost, ng, dc)
            %% define named indices into data matrices
            [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

            ipwl = find(cost(:, MODEL) == PW_LINEAR);   %% piece-wise linear costs
            ny = size(ipwl, 1);   %% number of piece-wise linear cost vars
            if dc
                nq = 0;    %% number of Qg variables
                q1 = [];
            else
                nq = ng;    %% number of Qg variables
                q1 = 1+ng;
            end
            [Ay, by] = makeAy(baseMVA, ng, cost, 1, q1, 1+ng+nq);
            p = struct('n', ny, 'i', ipwl, 'A', Ay, 'b', by);
        end

        function p = gen_cost_poly_params(obj, baseMVA, cost)
            %% define named indices into data matrices
            [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
            
            ng = size(cost, 1);
            have_quad_cost = 0;
            kg = []; cg = []; Qg = [];
            i0 = find(cost(:, MODEL) == POLYNOMIAL & cost(:, NCOST) == 1);   %% constant
            i1 = find(cost(:, MODEL) == POLYNOMIAL & cost(:, NCOST) == 2);   %% linear
            i2 = find(cost(:, MODEL) == POLYNOMIAL & cost(:, NCOST) == 3);   %% quadratic
            i3 = find(cost(:, MODEL) == POLYNOMIAL & cost(:, NCOST) > 3);    %% cubic or greater
            if ~isempty(i2) || ~isempty(i1) || ~isempty(i0)
                have_quad_cost = 1;
                kg = zeros(ng, 1);
                cg = zeros(ng, 1);
                if ~isempty(i2)
                    Qg = zeros(ng, 1);
                    Qg(i2) = 2 * cost(i2, COST) * baseMVA^2;
                    cg(i2) = cg(i2) + cost(i2, COST+1) * baseMVA;
                    kg(i2) = kg(i2) + cost(i2, COST+2);
                end
                if ~isempty(i1)
                    cg(i1) = cg(i1) + cost(i1, COST) * baseMVA;
                    kg(i1) = kg(i1) + cost(i1, COST+1);
                end
                if ~isempty(i0)
                    kg(i0) = kg(i0) + cost(i0, COST);
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

        function [A, l, u] = disp_load_constant_pf_constraint(obj, dm);
            %%-----  HACK ALERT  -----
            %% create a mpc with only online gens to call makeAvl
            mpc = dm.mpc;
            mpc.gen = mpc.gen(obj.on, :);
            %%-----  end of HACK  -----

            %% this should give incorrect results if mpc is not in 
            [A, l, u] = makeAvl(mpc);
        end
    end     %% methods
end         %% classdef
