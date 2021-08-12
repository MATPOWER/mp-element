classdef dme_load < dm_element
%DME_LOAD  MATPOWER data model class for load data

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus     %% bus index vector (all loads)
        Pd      %% active power demand (p.u.) for constant power loads that are on
        Qd      %% reactive power demand (p.u.) for constant power loads that are on
        Pd_i    %% active power demand (p.u.) for constant current loads that are on
        Qd_i    %% reactive power demand (p.u.) for constant current loads that are on
        Pd_z    %% active power demand (p.u.) for constant impedance loads that are on
        Qd_z    %% reactive power demand (p.u.) for constant impedance loads that are on
    end     %% properties

    methods
        %% constructor
        function obj = dme_load()
            obj@dm_element();   %% call parent constructor
            obj.name = 'load';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus', 'pd', 'qd', 'pd_i', 'qd_i', 'pd_z', 'qd_z', ...
                'p', 'q'});
        end

        function vars = export_vars(obj, task)
            switch task
                case 'PF'
                    vars = {};
                case 'CPF'
                    vars = {'pd', 'qd', 'pd_i', 'qd_i', 'pd_z', 'qd_z'};
                case 'OPF'
                    vars = {};
                otherwise
                    vars = 'all';
            end
        end

        function nr = count(obj, dm)
            nr = count@dm_element(obj, dm);
            obj.bus = obj.tab.source_uid;
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
            obj.Pd   = obj.tab.pd(obj.on) / dm.base_mva;
            obj.Qd   = obj.tab.qd(obj.on) / dm.base_mva;
            obj.Pd_i = obj.tab.pd_i(obj.on) / dm.base_mva;
            obj.Qd_i = obj.tab.qd_i(obj.on) / dm.base_mva;
            obj.Pd_z = obj.tab.pd_z(obj.on) / dm.base_mva;
            obj.Qd_z = obj.tab.qd_z(obj.on) / dm.base_mva;
        end

        function dm = parameterized(obj, dm, dmb, dmt, lam)
            load = dm.elements.load;
            b = dmb.elements.load.tab;      %% base load table
            t = dmt.elements.load.tab;      %% target load table

            load.tab.pd = b.pd + lam * (t.pd - b.pd);
            load.tab.qd = b.qd + lam * (t.qd - b.qd);
            load.tab.pd_i = b.pd_i + lam * (t.pd_i - b.pd_i);
            load.tab.qd_i = b.qd_i + lam * (t.qd_i - b.qd_i);
            load.tab.pd_z = b.pd_z + lam * (t.pd_z - b.pd_z);
            load.tab.qd_z = b.qd_z + lam * (t.qd_z - b.qd_z);
        end
    end     %% methods
end         %% classdef
