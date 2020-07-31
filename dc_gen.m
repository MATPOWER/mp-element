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
        function obj = add_zvars(obj, nm, mpc, idx)
            %% define constants
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
                MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
                QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

            ng = obj.nk;
            Pg   = mpc.gen(:, PG) / mpc.baseMVA;
            Pmin = mpc.gen(:, PMIN) / mpc.baseMVA;
            Pmax = mpc.gen(:, PMAX) / mpc.baseMVA;
            nm.add_var('z', 'Pg', ng, Pg, Pmin, Pmax);
        end

        function obj = build_params(obj, nm, mpc)
            build_params@mp_gen(obj, nm, mpc);      %% call parent
            ng = obj.nk;
            obj.K = -speye(ng);
        end

        function add_opf_constraints(obj, nm, om, mpc, mpopt)
            %% piecewise linear costs
            if obj.cost_pwl.ny
                om.add_lin_constraint('ycon', obj.cost_pwl.Ay, [], obj.cost_pwl.by, {'Pg', 'y'});
            end

            %% call parent
            add_opf_constraints@mp_gen(obj, nm, om, mpc, mpopt);
        end
    end     %% methods
end         %% classdef
