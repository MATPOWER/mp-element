classdef acc_bus < mp_bus & acc_model

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
        function obj = add_vvars(obj, asm, mpc, idx)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            nb = obj.nk;
            V0 = mpc.bus(:, VM) .* exp(1j * mpc.bus(:, VA) * pi/180);
            Vclim = 1.1 * mpc.bus(:, VMAX);
            asm.add_var('vr', 'Vr', nb, real(V0), -Vclim, Vclim);
            asm.add_var('vi', 'Vi', nb, imag(V0), -Vclim, Vclim);
        end
    end     %% methods
end         %% classdef
