classdef acsp_bus < mp_bus & acsp_model

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
        %% constructor
        function obj = acsp_bus(varargin)
            obj@mp_bus(varargin{:});
        end

        function obj = add_vvars(obj, asm, mpc, idx)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            nb = obj.nk;
            Va0   = mpc.bus(:, VA) * pi/180;
            Vm0   = mpc.bus(:, VM);
            Vamax = Inf(nb, 1);
            Vamin = -Inf(nb, 1);
            k = find(mpc.bus(:, BUS_TYPE) == REF);
            Vamax(k) = Va0(k);
            Vamin(k) = Va0(k);
            Vmin = mpc.bus(:, VMIN);
            Vmax = mpc.bus(:, VMAX);
            asm.add_var('va', 'Va', nb, Va0, Vamin, Vamax);
            asm.add_var('vm', 'Vm', nb, Vm0, Vmin, Vmax);
        end
    end     %% methods
end         %% classdef
