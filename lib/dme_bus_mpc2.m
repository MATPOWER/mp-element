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
            nb = obj.n;
            gbus = obj.i2on(gen_dme.bus(gen_dme.on));   %% buses of online gens
            bus = obj.get_table(dm);
            gen = gen_dme.get_table(dm);

            %% update bus types based on connected generator status
            %% gen connection matrix, element i, j is 1 if gen j @ bus i is ON
            Cg = sparse(gbus, (1:ng)', 1, nb, ng);
            bus_gen_status = Cg * ones(ng, 1);  %% num of gens ON at each bus
%             obj.isref = obj.isref & bus_gen_status;
              % above line would affect OPF (not just PF, CPF) where REF is
              % used only as angle reference and does not require an online gen
            obj.ispv = obj.ispv &  bus_gen_status;
            obj.ispq = obj.ispq | ~bus_gen_status;
%             obj.ensure_ref_bus();   %% pick a new ref bus if one does not exist

            %% initialize voltage from bus table
            obj.Va0 = bus(obj.on, VA) * pi/180;
            obj.Vm0 = bus(obj.on, VM);

            %% pull PV bus voltage magnitudes from mpc.gen(:, VG)
            vcb = ones(nb, 1);      %% create mask of voltage-controlled buses
            vcb(obj.ispq) = 0;      %% exclude PQ buses
            %% find indices of online gens at online v-c buses
            k = find(vcb(gbus));
            obj.Vm0(gbus(k)) = gen(gen_dme.on(k), VG);

            obj.Vmin = bus(obj.on, VMIN);
            obj.Vmax = bus(obj.on, VMAX);
        end

        function obj = set_bus_type_ref(obj, dm, idx)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;
            dm.mpc.bus(idx, BUS_TYPE) = REF;
            set_bus_type_ref@dme_bus(obj, dm, idx); %% call parent
        end

        function obj = set_bus_type_pv(obj, dm, idx)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;
            dm.mpc.bus(idx, BUS_TYPE) = PV;
            set_bus_type_pv@dme_bus(obj, dm, idx);  %% call parent
        end

        function obj = set_bus_type_pq(obj, dm, idx)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;
            dm.mpc.bus(idx, BUS_TYPE) = PQ;
            set_bus_type_pq@dme_bus(obj, dm, idx);  %% call parent
        end

        function obj = update(obj, dm, varargin)
            %% obj.update(dm, name1, val1, name2, val2, ...)
            %% obj.update(dm, idx, name1, val1, name2, val2, ...)

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            n = length(varargin);
            if rem(n, 2)    %% odd
                idx = obj.on(varargin{1});
                s = 2;      %% starting arg index
            else            %% even
                idx = obj.on;
                s = 1;      %% starting arg index
            end
            for k = s:2:n-1
                val = varargin{k+1};
                switch varargin{k}
                    case 'Va'
                        dm.mpc.bus(idx, VA) = val * 180/pi;
                    case 'Vm'
                        dm.mpc.bus(idx, VM) = val;
                    case 'bus_type'
                        obj.isref = (val == REF);   %% bus is ref?
                        obj.ispv  = (val == PV);    %% bus is PV?
                        obj.ispq  = (val == PQ);    %% bus is PQ?
                        dm.mpc.bus(idx, BUS_TYPE) = val;
                    case 'lamP'
                        dm.mpc.bus(idx, LAM_P) = val / dm.baseMVA;
                    case 'lamQ'
                        dm.mpc.bus(idx, LAM_Q) = val / dm.baseMVA;
                    case 'muVmin'
                        dm.mpc.bus(idx, MU_VMIN) = val;
                    case 'muVmax'
                        dm.mpc.bus(idx, MU_VMAX) = val;
                end
            end
        end
    end     %% methods
end         %% classdef
