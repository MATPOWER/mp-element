classdef dme_bus_ld_node_test < dme_bus_nld_node_test
%DME_BUS_LD_NODE_TEST  MATPOWER data model bus table for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020-2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        Pd      %% active power demand (p.u.) for loads at buses that are on
        Qd      %% reactive power demand (p.u.) for loads at buses that are on
        Gs      %% shunt conductance (p.u. active power demanded at
                %% V = 1.0 p.u.) for shunts at buses that are on
        Bs      %% shunt susceptance (p.u. reactive power injected at
                %% V = 1.0 p.u.) for shunts at buses that are on
    end     %% properties

    methods
        %% constructor
        function obj = dme_bus_ld_node_test()
            obj@dme_bus_nld_node_test();   %% call parent constructor
            obj.bus_eti = 2;        %% bus element type index, 2 => bus_ld
            obj.name = 'bus_ld';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( ...
                table_var_names@dme_bus_nld_node_test(obj), ...
                {'pd', 'qd', 'gs', 'bs'} );
        end

        function obj = build_params(obj, dm)
            obj = build_params@dme_bus(obj, dm);   %% call parent

            obj.Pd = obj.tab.pd(obj.on) / dm.baseMVA;
            obj.Qd = obj.tab.qd(obj.on) / dm.baseMVA;
            obj.Gs = obj.tab.gs(obj.on) / dm.baseMVA;
            obj.Bs = obj.tab.bs(obj.on) / dm.baseMVA;
        end
    end     %% methods
end         %% classdef
