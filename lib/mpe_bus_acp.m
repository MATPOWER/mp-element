classdef mpe_bus_acp < mpe_bus & mp_model_acp

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'bus';
%     end
    
    methods
        function obj = add_vvars(obj, nm, dm, idx)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
               MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
               QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            mpc = dm.mpc;
            nb = obj.nk;
            Va0   = mpc.bus(:, VA) * pi/180;
            Vm0   = mpc.bus(:, VM);

            %% pull gen bus voltages from mpc.gen(:, VG)
            on = find(mpc.gen(:, GEN_STATUS) > 0);  %% which generators are on?
            gbus = mpc.gen(on, GEN_BUS);            %% what buses are they at?
            vcb = ones(nb, 1);      %% create mask of voltage-controlled buses
            vcb(mpc.bus(:, BUS_TYPE) == PQ) = 0;    %% exclude PQ buses
            k = find(vcb(gbus));    %% v-c buses w/in-service gens
            Vm0(gbus(k)) = mpc.gen(on(k), VG);

            Vamin = -Inf(nb, 1);
            Vamax =  Inf(nb, 1);
            k = find(mpc.bus(:, BUS_TYPE) == REF);
            Vamin(k) = Va0(k);
            Vamax(k) = Va0(k);
            Vmin = mpc.bus(:, VMIN);
            Vmax = mpc.bus(:, VMAX);
            nm.add_var('va', 'Va', nb, Va0, Vamin, Vamax);
            nm.add_var('vm', 'Vm', nb, Vm0, Vmin, Vmax);
        end
    end     %% methods
end         %% classdef
