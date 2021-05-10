classdef dme_bus_ld_mpc2_node_test < dme_bus_mpc2
%DME_BUS_LD_MPC2_NODE_TEST  MATPOWER data model bus table for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
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
        function obj = dme_bus_ld_mpc2_node_test()
            obj@dme_bus_mpc2();     %% call parent constructor
            obj.name = 'bus_ld';
        end

        function obj = build_params(obj, dm)
            obj = build_params@dme_bus_mpc2(obj, dm);   %% call parent

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS] = idx_bus;

            tab = obj.get_table(dm);
            obj.Pd = tab(obj.on, PD) / dm.baseMVA;
            obj.Qd = tab(obj.on, QD) / dm.baseMVA;
            obj.Gs = tab(obj.on, GS) / dm.baseMVA;
            obj.Bs = tab(obj.on, BS) / dm.baseMVA;
        end
    end     %% methods
end         %% classdef
