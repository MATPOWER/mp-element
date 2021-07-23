classdef dme_gen < dm_element
%DME_GEN  Abstract base class for MATPOWER data model gen table

%   MATPOWER
%   Copyright (c) 1996-2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus     %% bus index vector (all gens)
        Pg0     %% initial active power (p.u.) for gens that are on
        Qg0     %% initial reactive power (p.u.) for gens that are on
        Vg      %% generator voltage setpoint
        Pmin    %% active power lower bound (p.u.) for gens that are on
        Pmax    %% active power upper bound (p.u.) for gens that are on
        Qmin    %% reactive power lower bound (p.u.) for gens that are on
        Qmax    %% reactive power upper bound (p.u.) for gens that are on
        cost_pg %% table for cost parameters for active power generation
        cost_qg %% table for cost parameters for reactive power generation
        pwl1    %% indices of single-block piecewise linear costs
                %% (automatically converted to linear cost)
    end     %% properties

    methods
        %% constructor
        function obj = dme_gen()
            obj@dm_element();   %% call parent constructor
            obj.name = 'gen';
            obj.table = 'gen';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus', 'vm_setpoint', 'pg_lb', 'pg_ub', 'qg_lb', 'qg_ub', ...
                'pg', 'qg', 'in_service', ...
                'startup_cost_cold', ...
                'pc1', 'pc2', 'qc1_lb', 'qc1_ub', 'qc2_lb', 'qc2_ub', ...
                'mu_pg_lb', 'mu_pg_ub', 'mu_qg_lb', 'mu_qg_ub'});
        end

        function vars = export_vars(obj, task)
            switch task
                case 'PF'
                    vars = {'pg', 'qg'};
                case 'CPF'
                    vars = {'pg', 'qg'};
                case 'OPF'
                    vars = {'vm_setpoint', 'pg', 'qg', 'mu_pg_lb', 'mu_pg_ub', 'mu_qg_lb', 'mu_qg_ub'};
                otherwise
                    vars = 'all';
            end
        end

        function obj = initialize(obj, dm)
            initialize@dm_element(obj, dm);    %% call parent

            %% get bus mapping info
            b2i = dm.elements.bus.ID2i;     %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            obj.bus = b2i(obj.tab.bus);
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elements.bus.status;    %% bus status

            %% update status of gens at isolated/offline buses
            obj.status = obj.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            baseMVA = dm.baseMVA;

            gen = obj.tab;

            %% get generator parameters
            obj.Pg0  = gen.pg(obj.on) / baseMVA;
            obj.Pmin = gen.pg_lb(obj.on) / baseMVA;
            obj.Pmax = gen.pg_ub(obj.on) / baseMVA;
            obj.Qg0  = gen.qg(obj.on) / baseMVA;
            obj.Qmin = gen.qg_lb(obj.on) / baseMVA;
            obj.Qmax = gen.qg_ub(obj.on) / baseMVA;
            obj.Vg   = gen.vm_setpoint(obj.on);
        end

        function obj = update(obj, dm, varargin)
            %% obj.update(dm, name1, val1, name2, val2, ...)
            %% obj.update(dm, idx, name1, val1, name2, val2, ...)

            %% define named indices into data matrices
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
            baseMVA = dm.baseMVA;

            n = length(varargin);
            if rem(n, 2)    %% odd
                idx = obj.on(varargin{1});
                s = 2;      %% starting arg index
            else            %% even
                idx = obj.on;
                s = 1;      %% starting arg index
            end
            for k = s:2:n-1
                val = varargin{k+1};
                switch varargin{k}
                    case 'Sg'
                        obj.tab.pg(idx) = real(val) * baseMVA;
                        obj.tab.qg(idx) = imag(val) * baseMVA;
                    case 'Pg'
                        obj.tab.pg(idx) = val * baseMVA;
                    case 'Qg'
                        obj.tab.qg(idx) = val * baseMVA;
                    case 'Vg'
                        obj.tab.vm_setpoint(idx) = val;
                    case 'muPmin'
                        obj.tab.mu_pg_lb(idx) = val / baseMVA;
                    case 'muPmax'
                        obj.tab.mu_pg_ub(idx) = val / baseMVA;
                    case 'muQmin'
                        obj.tab.mu_qg_lb(idx) = val / baseMVA;
                    case 'muQmax'
                        obj.tab.mu_qg_ub(idx) = val / baseMVA;
                end
            end
        end

        function [mn, mx, both] = violated_q_lims(obj, dm, mpopt)
            %% [mn, mx, both] = obj.violated_q_lims(dm, mpopt)
            %%  indices of online gens with violated Q lims

            gen = obj.tab;
            on = obj.on;

            %% find gens with violated Q constraints
            mx = find( gen.qg(on) > gen.qg_ub(on) + mpopt.opf.violation );
            mn = find( gen.qg(on) < gen.qg_lb(on) - mpopt.opf.violation );
            both = union(mx', mn')';    %% transposes handle fact that
                                        %% union of scalars is a row vector

            if ~isempty(both)   %% we have some Q limit violations
                %% first check for INFEASIBILITY
                %% find available online gens at REF and PV buses
                bus_dme = dm.elements.bus;
                gbus = bus_dme.i2on(obj.bus(obj.on));   %% buses of online gens
                remaining = find( bus_dme.isref(gbus) | bus_dme.ispv( gbus) );

                if length(both) == length(remaining) && ...
                        all(both == remaining) && (isempty(mx) || isempty(mn))
                    %% all remaining PV/REF gens are violating AND all are
                    %% violating same limit (all violating Qmin or all Qmax)
                    mn = [];
                    mx = [];
                else
                    %% one at a time?
                    if mpopt.pf.enforce_q_lims == 2
                        %% fix largest violation, ignore the rest
                        [junk, k] = max([gen.qg(mx) - gen.qg_ub(mx);
                                         gen.qg_lb(mn) - gen.qg(mn)]);
                        if k > length(mx)
                            mn = mn(k-length(mx));
                            mx = [];
                        else
                            mx = mx(k);
                            mn = [];
                        end
                    end
                end
            end
        end

        %%-----  OPF methods  -----
        function cost = opf_build_gen_cost_params(obj, dm, dc)
            baseMVA = dm.baseMVA;

            poly_p = obj.gen_cost_poly_params(baseMVA, obj.cost_pg);
            if dc || isempty(obj.cost_qg)
                pwl = obj.gen_cost_pwl_params(baseMVA, obj.cost_pg, obj.n, dc);
                poly_q = [];
            else
                poly_q = obj.gen_cost_poly_params(baseMVA, obj.cost_qg);

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
                pwl = obj.gen_cost_pwl_params(baseMVA, cost, obj.n, dc);
            end

            cost = struct( ...
                    'poly_p',   poly_p, ...
                    'poly_q',   poly_q, ...
                    'pwl',      pwl ...
                );
        end

        function p = gen_cost_pwl_params(obj, baseMVA, cost, ng, dc)
            ipwl = find(cost.pwl_n);    %% piece-wise linear costs
            ny = size(ipwl, 1);     %% number of piece-wise linear cost vars
            if dc
                nq = 0;     %% number of Qg variables
                q1 = [];
            else
                nq = ng;    %% number of Qg variables
                q1 = 1+ng;
            end

            %% from makeAy()
            ybas = 1+ng+nq;
            if ny == 0
                Ay = sparse([], [], [], 0, ybas+ny-1, 0);
                by = [];
            else
                %% if p(i),p(i+1),c(i),c(i+1) define one of the cost segments,
                %% then the corresponding constraint on Pg (or Qg) and Y is
                %%                                             c(i+1) - c(i)
                %%  Y   >=   c(i) + m * (Pg - p(i)),      m = ---------------
                %%                                             p(i+1) - p(i)
                %%
                %% this becomes   m * Pg - Y   <=   m*p(i) - c(i)

                %% form constraint matrix
                m = sum(cost.pwl_n(ipwl));  %% total number of cost points
                Ay = sparse([], [], [], m-ny, ybas+ny-1, 2*(m-ny)); 
                by = [];
                k = 1;
                j = 1;
                for i=ipwl'
                    ns = cost.pwl_n(i); %% # of cost points; segments = ns-1
                    p = cost.pwl_qty(i, 1:ns) / baseMVA;
                    c = cost.pwl_cost(i, 1:ns);
                    m = diff(c) ./ diff(p);         %% slopes for Pg (or Qg)
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

        function p = gen_cost_poly_params(obj, baseMVA, cost)
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
                    Qg(i2) = 2 * cost.poly_coef(i2, 3) * baseMVA^2;
                    cg(i2) = cg(i2) + cost.poly_coef(i2, 2) * baseMVA;
                    kg(i2) = kg(i2) + cost.poly_coef(i2, 1);
                end
                if ~isempty(i1)
                    cg(i1) = cg(i1) + cost.poly_coef(i1, 2) * baseMVA;
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

        %% from makeApq()
        function [Ah, uh, Al, ul, data] = pq_capability_constraint(obj, baseMVA);
            gen = obj.tab(obj.on, :);

            %% data dimensions
            ng = size(gen, 1);      %% number of dispatchable injections

            %% which generators require additional linear constraints
            %% (in addition to simple box constraints) on (Pg,Qg) to correctly
            %% model their PQ capability curves
            ipqh = find( obj.has_pq_cap(gen, 'U') );
            ipql = find( obj.has_pq_cap(gen, 'L') );
            npqh = size(ipqh, 1);   %% number of general PQ capability curves (upper)
            npql = size(ipql, 1);   %% number of general PQ capability curves (lower)

            %% make Ah if there is a need to add general PQ capability curves;
            %% use normalized coefficient rows so multipliers have right scaling
            %% in $$/pu
            if npqh > 0
                data.h = [gen.qc1_ub(ipqh)-gen.qc2_ub(ipqh), gen.pc2(ipqh)-gen.pc1(ipqh)];
                uh = data.h(:, 1) .* gen.pc1(ipqh) + data.h(:, 2) .* gen.qc1_ub(ipqh);
                for i=1:npqh,
                    tmp = norm(data.h(i,:));
                    data.h(i,:) = data.h(i, :) / tmp;
                    uh(i) = uh(i) / tmp;
                end
                Ah = sparse([1:npqh, 1:npqh]', [ipqh; ipqh+ng], ...
                            data.h(:), npqh, 2*ng);
                uh = uh / baseMVA;
            else
                data.h = [];
                Ah  = sparse(0, 2*ng);
                uh = [];
            end

            %% similarly Al
            if npql > 0
                data.l = [gen.qc2_lb(ipql)-gen.qc1_lb(ipql), gen.pc1(ipql)-gen.pc2(ipql)];
                ul= data.l(:, 1) .* gen.pc1(ipql) + data.l(:, 2) .* gen.qc1_lb(ipql) ;
                for i=1:npql,
                    tmp = norm(data.l(i, : ));
                    data.l(i, :) = data.l(i, :) / tmp;
                    ul(i) = ul(i) / tmp;
                end
                Al = sparse([1:npql, 1:npql]', [ipql; ipql+ng], ...
                            data.l(:), npql, 2*ng);
                ul = ul / baseMVA;
            else
                data.l = [];
                Al  = sparse(0, 2*ng);
                ul = [];
            end

            data.ipql = ipql;
            data.ipqh = ipqh;
        end

        %% from hasPQcap()
        function TorF = has_pq_cap(obj, gen, upper_lower)
            %% default value
            if nargin < 3
                upper_lower = 'B';  %% look at both top and bottom by default
            end

            %% for which gens is it specified
            k = find( gen.pc1 | gen.pc2 );
            ng = size(gen, 1);

            if isempty(k)
                TorF = zeros(ng, 1);
            else
                %% eliminate cases where QMIN = QMAX = QC
                kk = find(  gen.qg_lb(k) == gen.qg_ub(k) & ...
                            gen.qg_lb(k) == gen.qc1_ub(k) & ...
                            gen.qg_lb(k) == gen.qc1_lb(k) & ...
                            gen.qg_lb(k) == gen.qc2_ub(k) & ...
                            gen.qg_lb(k) == gen.qc2_lb(k) );
                k(kk) = [];

                %% check for errors in capability curve data
                if any( gen.pc1(k) >= gen.pc2(k) )
                    error('dme_gen/has_pq_cap: must have pc1 < pc2');
                end
                if any( gen.qc2_ub(k) <= gen.qc2_lb(k) & gen.qc1_ub(k) <= gen.qc1_lb(k) )
                    error('dme_gen/has_pq_cap: capability curve defines an empty set');
                end

                %% for which gens is it specified
                k = find( gen.pc1 ~= gen.pc2 );
                L = zeros(ng, 1);
                U = zeros(ng, 1);
                dPc = gen.pc2(k) - gen.pc1(k);

                if ~strcmp(upper_lower, 'U')    %% include lower constraint
                    dQc = gen.qc2_lb(k) - gen.qc1_lb(k);
                    Qmin_at_Pmin = gen.qc1_lb(k) + (gen.pg_lb(k) - gen.pc1(k)) .* ...
                        dQc ./ dPc;
                    Qmin_at_Pmax = gen.qc1_lb(k) + (gen.pg_ub(k) - gen.pc1(k)) .* ...
                        dQc ./ dPc;
                    L(k) = Qmin_at_Pmin > gen.qg_lb(k) | Qmin_at_Pmax > gen.qg_lb(k);
                end

                if ~strcmp(upper_lower, 'L')    %% include upper constraint
                    dQc = gen.qc2_ub(k) - gen.qc1_ub(k);
                    Qmax_at_Pmin = gen.qc1_ub(k) + (gen.pg_lb(k) - gen.pc1(k)) .* ...
                        dQc ./ dPc;
                    Qmax_at_Pmax = gen.qc1_ub(k) + (gen.pg_ub(k) - gen.pc1(k)) .* ...
                        dQc ./ dPc;
                    U(k) = Qmax_at_Pmin < gen.qg_ub(k) | Qmax_at_Pmax < gen.qg_ub(k);
                end

                TorF = L | U;
            end
        end

        %% from makeAvl()
        function [A, l, u] = disp_load_constant_pf_constraint(obj, dm);
            %% once this is working, let's eliminate this copying of the table
            baseMVA = dm.baseMVA;
            gen = obj.tab(obj.on, :);

            %% data dimensions
            ng = size(gen, 1);      %% number of dispatchable injections
            Pg = obj.Pg0;
            Qg = obj.Qg0;
            Pmin = obj.Pmin;
            Qmin = obj.Qmin;
            Qmax = obj.Qmax;

            ivl = find( obj.isload(obj.on) & (Qmin ~= 0 | Qmax ~= 0) );
            nvl  = size(ivl, 1);  %% number of dispatchable loads

            %% at least one of the Q limits must be zero (corresponding to Pmax == 0)
            if any( Qmin(ivl) ~= 0 & Qmax(ivl) ~= 0 )
                k = find(Qmin(ivl) ~= 0 & Qmax(ivl) ~= 0);
                gidx = obj.tab.uid(obj.on(ivl(k)));
                s = sprintf('Invalid qg limits for dispatchable load in row %d of gen table\n', gidx);
                error('dme_gen/disp_load_constant_pf_constraint: Either qg_lb or qg_ub must be equal to zero for each dispatchable load.\n%s', s);
            end

            %% Initial values of PG and QG must be consistent with specified
            %% power factor
            Qlim = (Qmin(ivl) == 0) .* Qmax(ivl) + ...
                (Qmax(ivl) == 0) .* Qmin(ivl);
            if any( abs( Qg(ivl) - Pg(ivl) .* Qlim ./ Pmin(ivl) ) > 1e-6 )
                k = find(abs( Qg(ivl) - Pg(ivl) .* Qlim ./ Pmin(ivl) ) > 1e-6);
                gidx = obj.tab.uid(obj.on(ivl(k)));
                s = sprintf('qg for dispatchable load in row %d of gen table must be pg * %g\n', [gidx Qlim ./ Pmin(ivl)]');
                error('dme_gen/disp_load_constant_pf_constraint: %s\n         %s\n         %s\n         %s\n%s', ...
                    'For a dispatchable load, pg and qg must be consistent', ...
                    'with the power factor defined by pg_lb and the relevant', ...
                    '(non-zero) qg_lb or qg_ub limit.', ...
                    'Note: Setting pg = qg = 0 satisfies this condition.', s);
            end

            %% make A, l, u, for l <= A * [Pg; Qg] <= u
            if nvl > 0
                xx = Pmin(ivl);
                yy = Qlim;
                pftheta = atan2(yy, xx);
                pc = sin(pftheta);
                qc = -cos(pftheta);
                ii = [ (1:nvl)'; (1:nvl)' ];
                jj = [ ivl; ivl+ng ];
                A = sparse(ii, jj, [pc; qc], nvl, 2*ng);
                l = zeros(nvl, 1);
                u = l;
            else
                A = sparse(0, 2*ng);
                l = [];
                u = [];
            end
        end

        function TorF = isload(obj, idx)
            if nargin > 1
                TorF = obj.tab.pg_lb(idx) < 0 & obj.tab.pg_ub(idx) == 0;
            else
                TorF = obj.tab.pg_lb < 0 & obj.tab.pg_ub == 0;
            end
        end
    end     %% methods
end         %% classdef
