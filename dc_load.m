classdef dc_load < mp_load & dc_model

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'load';
%     end
    
    methods
        %% constructor
        function obj = dc_load(varargin)
            obj@mp_load(varargin{:});
        end

        function obj = build_params(obj, asm, mpc)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            build_params@mp_load(obj, asm, mpc);    %% call parent
            obj.p = mpc.bus(:, PD) / mpc.baseMVA;
        end
    end     %% methods
end         %% classdef
