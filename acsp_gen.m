classdef acsp_gen < mp_gen & acsp_model

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
        function obj = acsp_gen(varargin)
            obj@mp_gen(varargin{:});
        end

        function obj = add_zvars(obj, asm, mpc)
            %% define constants
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            ng = obj.nk;
            Pg   = mpc.gen(:, PG) / mpc.baseMVA;
            Pmin = mpc.gen(:, PMIN) / mpc.baseMVA;
            Pmax = mpc.gen(:, PMAX) / mpc.baseMVA;
            Qg   = mpc.gen(:, QG) / mpc.baseMVA;
            Qmin = mpc.gen(:, QMIN) / mpc.baseMVA;
            Qmax = mpc.gen(:, QMAX) / mpc.baseMVA;
            asm.add_var('zr', 'Pg', ng, Pg, Pmin, Pmax);
            asm.add_var('zi', 'Qg', ng, Qg, Qmin, Qmax);
        end

        function obj = build_params(obj, asm, mpc)
            build_params@mp_gen(obj, asm, mpc);     %% call parent
            ng = obj.nk;
            obj.N = -speye(ng);
        end
    end     %% methods
end         %% classdef
