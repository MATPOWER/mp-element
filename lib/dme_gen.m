classdef dme_gen < dm_element
%DME_GEN  MATPOWER data model class for gen data

%   MATPOWER
%   Copyright (c) 1996-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus         %% bus index vector (all gens)
        pg_start    %% initial active power (p.u.) for gens that are on
        qg_start    %% initial reactive power (p.u.) for gens that are on
        vm_setpoint %% generator voltage setpoint for gens that are on
        pg_lb       %% active power lower bound (p.u.) for gens that are on
        pg_ub       %% active power upper bound (p.u.) for gens that are on
        qg_lb       %% reactive power lower bound (p.u.) for gens that are on
        qg_ub       %% reactive power upper bound (p.u.) for gens that are on
    end     %% properties

    methods
        %% constructor
        function obj = dme_gen()
            obj@dm_element();   %% call parent constructor
            obj.name = 'gen';
            obj.cxn_type = 'bus';
            obj.cxn_idx_prop = 'bus';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus', 'vm_setpoint', 'pg_lb', 'pg_ub', 'qg_lb', 'qg_ub', ...
                'pg', 'qg', 'in_service', ...
                'startup_cost_cold', ...
                'pc1', 'pc2', 'qc1_lb', 'qc1_ub', 'qc2_lb', 'qc2_ub'} );
        end

        function vars = export_vars(obj, task)
            vars = {'pg', 'qg'};
        end

        function TorF = have_cost(obj)
            TorF = 0;
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
            base_mva = dm.base_mva;

            gen = obj.tab;

            %% get generator parameters
            obj.pg_start = gen.pg(obj.on) / base_mva;
            obj.pg_lb = gen.pg_lb(obj.on) / base_mva;
            obj.pg_ub = gen.pg_ub(obj.on) / base_mva;
            obj.qg_start  = gen.qg(obj.on) / base_mva;
            obj.qg_lb = gen.qg_lb(obj.on) / base_mva;
            obj.qg_ub = gen.qg_ub(obj.on) / base_mva;
            obj.vm_setpoint = gen.vm_setpoint(obj.on);
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
                %% bus types for buses with online gens
                bt = bus_dme.type(bus_dme.i2on(obj.bus(obj.on)));
                remaining = find( bt == NODE_TYPE.REF | bt == NODE_TYPE.PV );

                if length(both) == length(remaining) && ...
                        all(both == remaining) && (isempty(mx) || isempty(mn))
                    %% all remaining PV/REF gens are violating AND all are
                    %% violating same limit (all violating qg_lb or all qg_ub)
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

        function TorF = isload(obj, idx)
            if nargin > 1
                TorF = obj.tab.pg_lb(idx) < 0 & obj.tab.pg_ub(idx) == 0;
            else
                TorF = obj.tab.pg_lb < 0 & obj.tab.pg_ub == 0;
            end
        end
    end     %% methods
end         %% classdef
