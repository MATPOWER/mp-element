classdef dme_load_mpc2 < dme_load & dm_format_mpc2
%DME_LOAD_MPC2  MATPOWER data model load table for MATPOWER case format v2

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
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD] = idx_bus;

            %% get bus indices
            tab = obj.get_table(dm);
            obj.busidx = find(tab(:, PD) | tab(:, QD));

            %% number of loads
            nr = length(obj.busidx);
            obj.nr = nr;
        end

        function obj = build_params(obj, dm)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD] = idx_bus;
            baseMVA = dm.mpc.baseMVA;

            tab = obj.get_table(dm);
            obj.Pd = tab(obj.busidx, PD) / baseMVA;
            obj.Qd = tab(obj.busidx, QD) / baseMVA;
        end
    end     %% methods
end         %% classdef
