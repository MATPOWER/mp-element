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
            baseMVA = dm.mpc.baseMVA;

            tab = obj.get_table(dm);
            busidx = find(tab(:, GS) | tab(:, BS));
            obj.Gs = tab(busidx, GS) / baseMVA;
            obj.Bs = tab(busidx, BS) / baseMVA;
            %% temporarily store bus indices, until all indexing is
            %% finished and we can convert back to IDs
            obj.busID = busidx;
            nr = length(busidx);
            obj.nr = nr;
        end

        function obj = initialize(obj, dm)
            obj = initialize@dme_shunt(obj, dm);    %% call parent

            dme_bus = dm.elm_by_name('bus');
            obj.busID = dme_bus.ID(obj.busID);      %% convert bus idx to ID
        end
    end     %% methods
end         %% classdef
