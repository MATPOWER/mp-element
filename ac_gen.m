classdef ac_gen < mp_gen & acsp_model

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        cost_poly_q = struct();
    end
    
    methods
        %% constructor
        function obj = ac_gen(varargin)
            obj@mp_gen(varargin{:});
        end

        function obj = add_zvars(obj, asm, mpc, idx)
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

        function build_gen_cost_params(obj, mpc, mpopt, pcost, qcost)
            [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
            ng = obj.nk;
            if nargin < 5
                [pcost qcost] = pqcost(mpc.gencost, ng);
            end

            %% call parent (with pcost) to handle active power costs
            build_gen_cost_params@mp_gen(obj, mpc, mpopt, pcost);
            
            %% find/prepare polynomial generator costs for reactive power
            if ~isempty(qcost)
                iq0 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) == 1);   %% constant
                iq1 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) == 2);   %% linear
                iq2 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) == 3);   %% quadratic
                iq3 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) > 3);    %% cubic or greater
                if ~isempty(iq2) || ~isempty(iq1) || ~isempty(iq0)
                    kqg = zeros(ng, 1);
                    cqg = zeros(ng, 1);
                    if ~isempty(iq2)
                        Qqg = zeros(ng, 1);
                        Qqg(iq2) = 2 * qcost(iq2, COST) * mpc.baseMVA^2;
                        cqg(iq2) = cqg(iq2) + qcost(iq2, COST+1) * mpc.baseMVA;
                        kqg(iq2) = kqg(iq2) + qcost(iq2, COST+2);
                    else
                        Qqg = [];   %% no quadratic terms
                    end
                    if ~isempty(iq1)
                        cqg(iq1) = cqg(iq1) + qcost(iq1, COST) * mpc.baseMVA;
                        kqg(iq1) = kqg(iq1) + qcost(iq1, COST+1);
                    end
                    if ~isempty(iq0)
                        kqg(iq0) = kqg(iq0) + qcost(iq0, COST);
                    end
                end
                obj.cost_poly_q = struct( ...
                        'iq0', iq0, ...
                        'iq1', iq1, ...
                        'iq2', iq2, ...
                        'iq3', iq3, ...
                        'kqg', kqg, ...
                        'cqg', cqg, ...
                        'Qqg', Qqg ...
                    );
                if ~isempty(iq3)
                    error('mp_gen/build_gen_cost_params: polynomial generator costs greater than quadratic order not yet implemented');
                end
            end
        end
    end     %% methods
end         %% classdef
