classdef dc_bus < mp_bus & mp_model_dc

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'bus';
%     end
    
    methods
        function obj = add_vvars(obj, nm, mpc, idx)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            nb = obj.nk;
            Va0   = mpc.bus(:, VA) * pi/180;
            Vamax = Inf(nb, 1);
            k = find(mpc.bus(:, BUS_TYPE) == REF);
            Vamax(k) = Va0(k);
            nm.add_var('va', 'Va', nb, Va0, -Vamax, Vamax);
        end
    end     %% methods
end         %% classdef
