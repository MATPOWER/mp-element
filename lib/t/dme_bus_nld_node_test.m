classdef dme_bus_nld_node_test < dme_bus
%DME_BUS_NLD_NODE_TEST  MATPOWER data model class for non-load bus data for T_NODE_TEST

%   MATPOWER
%   Copyright (c) 2020-2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus_eti %% bus element type index, 1 = bus_nld, 2 = bus_ld
        bus    %% indices into bus matrix (all rows) for this type of bus
    end     %% properties

    methods
        %% constructor
        function obj = dme_bus_nld_node_test()
            obj@dme_bus();     %% call parent constructor
            obj.bus_eti = 1;        %% bus element type index, 1 => bus_nld
            obj.name = 'bus_nld';
        end

        function nr = count(obj, dm)
            nr = count@dme_bus(obj, dm);
            if nr
                obj.bus = obj.tab.source_uid;
            end
        end

        function [gbus, ig] = gbus_vector(obj, gen_dme)
            %% buses of online gens
            ig = find(gen_dme.bus_etv == obj.bus_eti);
            gbus = obj.i2on(gen_dme.bus(gen_dme.on(gen_dme.bus_etv == obj.bus_eti)));
        end
    end     %% methods
end         %% classdef
