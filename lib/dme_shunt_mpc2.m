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
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS] = idx_bus;

            %% get bus indices
            tab = obj.get_table(dm);
            obj.bus = find(tab(:, GS) | tab(:, BS));

            %% number of shunts
            nr = length(obj.bus);
            obj.nr = nr;
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elm_by_name('bus').status;  %% bus status

            %% update status of gens at isolated/offline buses
            obj.status = obj.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@dme_shunt(obj, dm);
        end

        function obj = build_params(obj, dm)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS] = idx_bus;
            baseMVA = dm.mpc.baseMVA;

            tab = obj.get_table(dm);
            obj.Gs = tab(obj.bus, GS) / baseMVA;
            obj.Bs = tab(obj.bus, BS) / baseMVA;
        end
    end     %% methods
end         %% classdef
