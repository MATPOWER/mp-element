classdef ac_load < mp_load & ac_model

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'load';
%     end

    methods
        function k = load_bus(obj, mpc)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            k = find(mpc.bus(:, PD) | mpc.bus(:, QD));
        end

        function obj = build_params(obj, asm, mpc)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            build_params@mp_load(obj, asm, mpc);    %% call parent

            obj.s = (mpc.bus(obj.busidx, PD) + ...
                1j * mpc.bus(obj.busidx, QD)) / mpc.baseMVA;    %% vector of complex power demand
        end
    end     %% methods
end         %% classdef
