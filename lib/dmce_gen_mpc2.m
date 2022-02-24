classdef dmce_gen_mpc2 < dmc_element_mpc2 % & dmce_gen
%DMCE_GEN_MPC2  Data model converter for gen elements for MATPOWER case v2.

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function name = name(obj)
            name = 'gen';
        end

        function table = table(obj)
            table = 'gen';
        end

        function vmap = table_var_map(obj, var_names, mpc)
            vmap = table_var_map@dmc_element_mpc2(obj, var_names, mpc);

            %% define named indices into data matrices
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            %% map type for each name (default mapping is -1)
            vmap.uid.type           = 3;    %% consecutive IDs, starting at 1
            vmap.name.type          = 2;    %% empty char
            vmap.source_uid.type    = 2;    %% empty char
            vmap.startup_cost_cold.type = 5;%% function call

            sci_fcn = @(ob, vn, nr, r, d)start_cost_import(ob, vn, nr, r, d);

            %% map arguments for each name
           %vmap.uid.args           = [];
           %vmap.name.args          = [];
            vmap.status.args        = GEN_STATUS;
           %vmap.source_uid.args    = [];
            vmap.bus.args           = GEN_BUS;
            vmap.vm_setpoint.args   = VG;
            vmap.pg_lb.args         = PMIN;
            vmap.pg_ub.args         = PMAX;
            vmap.qg_lb.args         = QMIN;
            vmap.qg_ub.args         = QMAX;
            vmap.pc1.args           = PC1;
            vmap.pc2.args           = PC2;
            vmap.qc1_lb.args        = QC1MIN;
            vmap.qc1_ub.args        = QC1MAX;
            vmap.qc2_lb.args        = QC2MIN;
            vmap.qc2_ub.args        = QC2MAX;
            vmap.pg.args            = PG;
            vmap.qg.args            = QG;
            vmap.in_service.args    = GEN_STATUS;
            vmap.startup_cost_cold.args = {sci_fcn};
            vmap.mu_pg_lb.args      = MU_PMIN;
            vmap.mu_pg_ub.args      = MU_PMAX;
            vmap.mu_qg_lb.args      = MU_QMIN;
            vmap.mu_qg_ub.args      = MU_QMAX;
        end

        function vals = start_cost_import(obj, vn, nr, r, mpc)
            %% define named indices into data matrices
            [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
            if isfield(mpc, 'gencost')
                if nr && isempty(r)
                    vals = mpc.gencost(1:nr, STARTUP);
                else
                    vals = mpc.gencost(r, STARTUP);
                end
            else
                vals = zeros(nr, 1);
            end
        end

        function tab = create_cost_table(obj, gencost);
            %% define named indices into data matrices
            [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

            if isempty(gencost)
                npoly = []; p = []; npwl = []; qty = []; cst = [];
            else
                ipwl = find(gencost(:, MODEL) == PW_LINEAR);
                ipol = find(gencost(:, MODEL) == POLYNOMIAL);
                nr = size(gencost, 1);

                %% get dimensions
                npoly = zeros(nr, 1);
                npwl = zeros(nr, 1);
                if isempty(ipol)
                    maxNpoly = 1;
                else
                    polycost = gencost(ipol, :);
                    npoly(ipol) = polycost(:, NCOST);
                    maxNpoly = max(npoly(ipol));
                    minNpoly = min(npoly(ipol));
                end
                if isempty(ipwl)
                    maxNpwl = 0;
                else
                    npwl(ipwl) = gencost(ipwl, NCOST);
                    maxNpwl = max(npwl);
                end

                %% initialize cost parameters
                p = zeros(nr, maxNpoly);
                qty = zeros(nr, maxNpwl);
                cst = zeros(nr, maxNpwl);

                %% polynomial costs, form coefficient matrix where 1st column
                %% is constant term, 2nd linear, etc.
                if ~isempty(ipol)
                    for n = minNpoly:maxNpoly
                        k = find(npoly(ipol) == n); %% cost with n coefficients
                        p(ipol(k), 1:n) = polycost(k, (COST+n-1):-1:COST);
                    end
                end

                %% piecewise linear costs
                if ~isempty(ipwl)
                    m = COST-1 + 2*maxNpwl;
                    qty(ipwl, :) = gencost(ipwl, COST:2:(m-1));
                    cst(ipwl, :) = gencost(ipwl, COST+1:2:m);
                end
            end
            var_names = {'poly_n', 'poly_coef', 'pwl_n', 'pwl_qty', 'pwl_cost'};
            if have_feature('table')
                tab = table(npoly, p, npwl, qty, cst, 'VariableNames', var_names);
            else
                tab = mp_table(npoly, p, npwl, qty, cst, 'VariableNames', var_names);
            end
        end

        function dme = import(obj, dme, mpc)
            %% call parent
            dme = import@dmc_element_mpc2(obj, dme, mpc);
            nr = size(dme.tab, 1);

            %% import gencost
            if dme.have_cost() && isfield(mpc, 'gencost') && nr
                %% define named indices into data matrices
                [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

                %% convert single-block piecewise-linear costs into linear polynomial cost
                gencost = mpc.gencost;
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
                    dme.pwl1 = pwl1;
                end
                [pcost, qcost] = pqcost(gencost, nr);
                dme.cost_pg = create_cost_table(obj, pcost);
                dme.cost_qg = create_cost_table(obj, qcost);
            end
        end

        function mpc = export(obj, dme, mpc, varargin)
            %% call parent
            mpc = export@dmc_element_mpc2(obj, dme, mpc, varargin{:});

            %% export gencost
            
        end
    end     %% methods
end         %% classdef
