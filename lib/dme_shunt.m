classdef dme_shunt < dm_element
%DME_SHUNT  MATPOWER data model class for shunt data

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus     %% bus index vector (all shunts)
        Gs      %% shunt conductance (p.u. active power demanded at
                %% V = 1.0 p.u.) for shunts that are on
        Bs      %% shunt susceptance (p.u. reactive power injected at
                %% V = 1.0 p.u.) for shunts that are on
    end     %% properties

    methods
        %% constructor
        function obj = dme_shunt()
            obj@dm_element();   %% call parent constructor
            obj.name = 'shunt';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus', 'gs', 'bs', 'p', 'q'});
        end

        function vars = export_vars(obj, task)
            switch task
                case 'PF'
                    vars = {};
                case 'CPF'
                    vars = {'gs', 'bs'};
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
            obj.Gs = obj.tab.gs(obj.on) / dm.baseMVA;
            obj.Bs = obj.tab.bs(obj.on) / dm.baseMVA;
        end

        function dm = parameterized(obj, dm, dmb, dmt, lam)
            shunt = dm.elements.shunt;
            b = dmb.elements.shunt.tab;     %% base shunt table
            t = dmt.elements.shunt.tab;     %% target shunt table

            shunt.tab.gs = b.gs + lam * (t.gs - b.gs);
            shunt.tab.bs = b.bs + lam * (t.bs - b.bs);
        end
    end     %% methods
end         %% classdef
