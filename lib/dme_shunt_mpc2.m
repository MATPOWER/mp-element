classdef dme_shunt_mpc2 < dme_shunt & dm_format_mpc2
%DME_SHUNT_MPC2  MATPOWER data model shunt table for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function nr = count(obj, dm)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS] = idx_bus;

            %% get bus indices
            tab = obj.get_table(dm);
            obj.busidx = find(tab(:, GS) | tab(:, BS));

            %% number of shunts
            nr = length(obj.busidx);
            obj.nr = nr;
        end

        function obj = build_params(obj, dm)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS] = idx_bus;
            baseMVA = dm.mpc.baseMVA;

            tab = obj.get_table(dm);
            obj.Gs = tab(obj.busidx, GS) / baseMVA;
            obj.Bs = tab(obj.busidx, BS) / baseMVA;
        end
    end     %% methods
end         %% classdef
