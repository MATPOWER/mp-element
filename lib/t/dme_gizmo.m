classdef dme_gizmo < dm_element
%DME_GIZMO  MATPOWER data model class for gizmo data

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus1        %% bus index vector for port 1
        bus2        %% bus index vector for port 2
        bus3        %% bus index vector for port 3
    end     %% properties

    methods
        %% constructor
        function obj = dme_gizmo()
            obj@dm_element();   %% call parent constructor
            obj.name = 'gizmo';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus_1', 'bus_2', 'bus_3', 'Y1r', 'Y1i', 'Y2r', 'Y2i', ...
                'Lr', 'Li', 'Ir', 'Ii', 'M1r', 'M1i', 'M2r', 'M2i', ...
                'Nr', 'Ni', 'Sr', 'Si', 'Zr1', 'Zi1', 'Zr2', 'Zi2'});
        end

        function obj = initialize(obj, dm)
            initialize@dm_element(obj, dm);     %% call parent

            %% get bus mapping info
            b2i = dm.elements.bus.ID2i;         %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            obj.bus1 = b2i(obj.tab.bus_1);
            obj.bus2 = b2i(obj.tab.bus_2);
            obj.bus3 = b2i(obj.tab.bus_3);
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elements.bus.status;        %% bus status

            %% update status of gizmoes connected to isolated/offline buses
            obj.status = obj.status & bs(obj.bus1) & ...
                                      bs(obj.bus2) & ...
                                      bs(obj.bus3);

            %% call parent to fill in on/off
            update_status@dm_element(obj, dm);
        end
    end     %% methods
end         %% classdef
