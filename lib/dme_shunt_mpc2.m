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

            tab = obj.get_table(dm);
            bidx = obj.bus(obj.on);             %% buses of shunts
            obj.Gs = tab(bidx, GS) / dm.baseMVA;
            obj.Bs = tab(bidx, BS) / dm.baseMVA;
        end

        function dm = parameterized(obj, dm, dmb, dmt, lam)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            cols = [GS BS];
            dm.mpc.bus(:, cols) = dmb.mpc.bus(:, cols) + ...
                lam * (dmt.mpc.bus(:, cols) - dmb.mpc.bus(:, cols));
        end
    end     %% methods
end         %% classdef
