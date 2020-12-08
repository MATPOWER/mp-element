classdef dme_bus_mpc2 < dme_bus & dm_format_mpc2
%DME_BUS  MATPOWER data model bus table for MATPOWER case format v2

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
        %% constructor
        function obj = dme_bus_mpc2()
            obj@dme_bus();      %% call parent constructor

            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;
            obj.id_col = BUS_I;
            obj.st_col = BUS_TYPE;
        end

        function status = get_status(obj, dm)
            %% overrides dm_element/get_status()
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;

            %% check that all buses have a valid BUS_TYPE
            tab = obj.get_table(dm);
            bt = tab(:, obj.st_col);
            err = find(~(bt == PQ | bt == PV | bt == REF | bt == NONE));
            if ~isempty(err)
                error('dme_bus_mpc2/get_status: bus %d has an invalid BUS_TYPE', err);
            end
            obj.isref = (bt == REF);    %% bus is ref?
            obj.ispv  = (bt == PV);     %% bus is PV?
            obj.ispq  = (bt == PQ);     %% bus is PQ?
            status = (bt ~= NONE);      %% bus status
            obj.status = status;
        end

        function obj = build_params(obj, dm)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
               MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
               QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            gen_dme = dm.elm_by_name('gen');
            bus = obj.get_table(dm);
            gen = gen_dme.get_table(dm);

            %% initialize voltage from bus table
            Va = bus(:, VA) * pi/180;
            Vm = bus(:, VM);

            %% pull PV bus voltage magnitudes from mpc.gen(:, VG)
            gbus = gen_dme.bus(gen_dme.on);     %% buses of online gens
            vcb = ones(obj.nr, 1);  %% create mask of voltage-controlled buses
            vcb(bus(:, BUS_TYPE) == PQ) = 0;    %% exclude PQ buses
            %% find indices of online at online v-c buses
            k = find(obj.status(gbus) & vcb(gbus));
            Vm(gbus(k)) = gen(gen_dme.on(k), VG);

            %% set initialize 
            obj.Vm0 = Vm(obj.on);
            obj.Va0 = Va(obj.on);
            obj.Vmin = bus(obj.on, VMIN);
            obj.Vmax = bus(obj.on, VMAX);
        end

        function btv = bus_types(obj, dm)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            gen_dme = dm.elm_by_name('gen');
            ng = gen_dme.n;

            %% generator connection matrix
            %% element i, j is 1 if, generator j at bus i is ON
            Cg = sparse(gen_dme.bus, (1:ng)', gen_dme.status, obj.n, ng);

            %% number of generators at each bus that are ON
            bus_gen_status = Cg * ones(ng, 1);

            isref = obj.isref & bus_gen_status;
            ispv  = obj.ispv  & bus_gen_status;
            ispq  = obj.ispq | ~bus_gen_status;

            %% pick a new reference bus if for some reason there is none
            %% (may have been shut down)
            if ~any(isref)
                k = find(ispv, 1);  %% find the first PV bus ...
                if isempty(k)
                    error('dme_bus_mpc2/bus_types: must have at least one REF or PV bus');
                end
                ispv(k)  = 0;       %% ...and change it to a REF bus
                isref(k) = 1;
            end

            %% package up bus type vector
            btv = isref * REF + ispv * PV + ispq * PQ;
        end
    end     %% methods
end         %% classdef
