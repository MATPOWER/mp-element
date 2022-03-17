classdef dme_buslink < dm_element
%DME_BUSLINK  MATPOWER data model class for 1-to-3-phase buslink data

%   MATPOWER
%   Copyright (c) 1996-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus         %% bus index vector (all buslinks)
        bus3p       %% bus3p index vector (all buslinks)
        pg1_start   %% initial phase 1 active power (p.u.) for buslinks that are on
        pg2_start   %% initial phase 2 active power (p.u.) for buslinks that are on
        pg3_start   %% initial phase 3 active power (p.u.) for buslinks that are on
        qg1_start   %% initial phase 1 reactive power (p.u.) for buslinks that are on
        qg2_start   %% initial phase 2 reactive power (p.u.) for buslinks that are on
        qg3_start   %% initial phase 3 reactive power (p.u.) for buslinks that are on
    end     %% properties

    methods
        function name = name(obj)
            name = 'buslink';
        end

        function label = label(obj)
            label = 'Bus Link';
        end

        function label = labels(obj)
            label = 'Bus Links';
        end

        function name = cxn_type(obj)
            name = {'bus', 'bus3p'};
        end

        function name = cxn_idx_prop(obj)
            name = {'bus', 'bus3p'};
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus', 'bus3p'});
        end

%         function vars = export_vars(obj, task)
%             vars = {};
%         end

        function obj = initialize(obj, dm)
            initialize@dm_element(obj, dm); %% call parent

            %% get bus mapping info
            b2i  = dm.elements.bus.ID2i;    %% bus num to idx mapping
            b2i3 = dm.elements.bus3p.ID2i;  %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            obj.bus   = b2i( obj.tab.bus);
            obj.bus3p = b2i3(obj.tab.bus3p);
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs  = dm.elements.bus.tab.status;   %% bus status
            bs3 = dm.elements.bus3p.tab.status; %% bus3p status

            %% update status of buslinks at isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.bus) & bs3(obj.bus3p);

            %% call parent to fill in on/off
            update_status@dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            %% check for matching base_kv
            base_kv1 = dm.elements.bus.tab.base_kv(obj.bus);
            base_kv3 = dm.elements.bus3p.tab.base_kv(obj.bus3p);
            if any(base_kv1 ~= base_kv3)
                error('dme_buslink/build_params: buslink objects can only link buses with identical base_kv');
            end

            %% check for at least one PQ bus per buslink
            type1 = dm.elements.bus.tab.type(obj.bus);
            type3 = dm.elements.bus3p.tab.type(obj.bus3p);
            if any(type1 ~= NODE_TYPE.PQ & type3 ~= NODE_TYPE.PQ)
                error('dme_buslink/build_params: at lease one of the buses linked by a buslink must be of type PQ');
            end

            %% check for balanced voltage magnitudes for REF and PV buses
            %% on 3-phase side
            refpv = (type3 == NODE_TYPE.REF | type3 == NODE_TYPE.PV);
            vm1_ref_pv = dm.elements.bus3p.tab.vm1(obj.bus3p(refpv));
            vm2_ref_pv = dm.elements.bus3p.tab.vm2(obj.bus3p(refpv));
            vm3_ref_pv = dm.elements.bus3p.tab.vm3(obj.bus3p(refpv));
            if any(vm1_ref_pv ~= vm2_ref_pv | vm2_ref_pv ~= vm3_ref_pv)
                error('dme_buslink/build_params: buslink objects can only link to REF or PV buses with balanced voltage magnitudes');
            end
        end
    end     %% methods
end         %% classdef
