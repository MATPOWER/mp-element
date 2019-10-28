classdef dc_gen < mp_gen & dc_model

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'gen';
%     end
    
    methods
        %% constructor
        function obj = dc_gen(varargin)
            obj@mp_gen(varargin{:});
        end

        function obj = add_zvars(obj, asm, mpc)
            %% define constants
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            ng = obj.nk;
            Pg   = mpc.gen(:, PG);
            Pmin = mpc.gen(:, PMIN);
            Pmax = mpc.gen(:, PMAX);
            asm.add_var('z', 'Pg', ng, Pg, Pmin, Pmax);
        end

        function obj = build_params(obj, asm, mpc)
            build_params@mp_gen(obj, asm, mpc);     %% call parent
            ng = obj.nk;
            obj.K = -speye(ng);
        end
    end     %% methods
end         %% classdef
