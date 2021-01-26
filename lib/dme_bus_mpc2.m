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

            %% define constants & named indices into data matrices
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

            %% temporarily set bus type properties with dimensions for all buses
            %% (reduced for online buses only in update_status())
            obj.isref = (bt == REF);    %% bus is ref?
            obj.ispv  = (bt == PV);     %% bus is PV?
            obj.ispq  = (bt == PQ);     %% bus is PQ?

            status = (bt ~= NONE);      %% bus status
            obj.status = status;
        end

        function obj = update_status(obj, dm)
            %% call parent to fill in on/off
            update_status@dme_bus(obj, dm);

            %% update bus type properties so they correspond
            %% to online buses only
            obj.isref = obj.isref(obj.on);
            obj.ispv  = obj.ispv(obj.on);
            obj.ispq  = obj.ispq(obj.on);
        end

        function obj = build_params(obj, dm)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
               MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
               QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            gen_dme = dm.elm_by_name('gen');
            ng = gen_dme.n;
            on = gen_dme.on;
            bus = obj.get_table(dm);
            gen = gen_dme.get_table(dm);

            %% update bus types based on connected generator status
            %% gen connection matrix, element i, j is 1 if gen j @ bus i is ON
            Cg = sparse(gen_dme.bus(on), (1:ng)', 1, obj.n, ng);
            bus_gen_status = Cg * ones(ng, 1);  %% num of gens ON at each bus
%             obj.isref = obj.isref & bus_gen_status;
              % above line would affect OPF (not just PF, CPF) where REF is
              % used only as angle reference and does not require an online gen
            obj.ispv  = obj.ispv  & bus_gen_status;
            obj.ispq  = obj.ispq | ~bus_gen_status;
%             obj.ensure_ref_bus();   %% pick a new ref bus if one does not exist

            %% initialize voltage from bus table
            Va = bus(:, VA) * pi/180;
            Vm = bus(:, VM);

            %% pull PV bus voltage magnitudes from mpc.gen(:, VG)
            gbus = gen_dme.bus(gen_dme.on);     %% buses of online gens
            vcb = ones(obj.nr, 1);  %% create mask of voltage-controlled buses
            vcb(obj.ispq) = 0;      %% exclude PQ buses
            %% find indices of online at online v-c buses
            k = find(obj.status(gbus) & vcb(gbus));
            Vm(gbus(k)) = gen(gen_dme.on(k), VG);

            obj.Vm0 = Vm(obj.on);
            obj.Va0 = Va(obj.on);
            obj.Vmin = bus(obj.on, VMIN);
            obj.Vmax = bus(obj.on, VMAX);
        end


        function obj = update(obj, dm, Va, Vm, lamP, lamQ, muVmin, muVmax)
            %% obj.update(dm, Va)
            %% obj.update(dm, Va, Vm)
            %% obj.update(dm, Va, [], lamP)
            %% obj.update(dm, Va, Vm, lamP, lamQ, muVmin, muVmax)

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            if nargin < 4 || isempty(Vm)
                Vm = 1;
            end

            dm.mpc.bus(obj.on, VA) = Va * 180/pi;
            dm.mpc.bus(obj.on, VM) = Vm;
            
            if nargin > 4
                dm.mpc.bus(obj.on, LAM_P) = lamP / dm.baseMVA;
                if nargin > 5
                    dm.mpc.bus(obj.on, LAM_Q) = lamQ / dm.baseMVA;
                    if nargin > 6
                        dm.mpc.bus(obj.on, MU_VMIN) = muVmin;
                        dm.mpc.bus(obj.on, MU_VMAX) = muVmax;
                    end
                end
            end
        end
    end     %% methods
end         %% classdef
