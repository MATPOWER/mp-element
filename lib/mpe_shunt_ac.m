classdef mpe_shunt_ac < mpe_shunt% & mp_model_ac

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'shunt';
%     end

    methods
        function obj = build_params(obj, nm, dm)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            build_params@mpe_shunt(obj, nm, dm);   %% call parent

            mpc = dm.mpc;
            nsh = obj.nk;
            Ysh = (mpc.bus(obj.busidx, GS) + ...
                1j * mpc.bus(obj.busidx, BS)) / mpc.baseMVA;    %% vector of shunt admittances
            obj.Y = sparse(1:nsh, 1:nsh, Ysh, nsh, nsh);
        end
    end     %% methods
end         %% classdef
