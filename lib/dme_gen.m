classdef dme_gen < dm_element
%DME_GEN  Abstract base class for MATPOWER data model gen table

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
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
        pcost   %% active power cost parameters for gens that are on
        qcost   %% reactive power cost parameters for gens that are on
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
                'mu_pg_lb', 'mu_pg_ub', 'mu_qg_lb', 'mu_qg_ub'});
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

            %% get gen cost parameters
            obj.build_cost_params(dm);
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
    end     %% methods
end         %% classdef
