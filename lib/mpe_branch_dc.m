classdef mpe_branch_dc < mpe_branch & mp_model_dc

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'branch';
%     end
    
    methods
        function obj = build_params(obj, nm, dm)
            build_params@mpe_branch(obj, nm, dm);   %% call parent

            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

            branch = dm.mpc.branch;
            nl = obj.nk;

            stat = branch(:, BR_STATUS);    %% ones at in-service branches
            tap = ones(nl, 1);              %% default tap ratio = 1
            i = find(branch(:, TAP));       %% indices of non-zero tap ratios
            tap(i) = branch(i, TAP);        %% assign non-zero tap ratios

            b = stat ./ branch(:, BR_X);    %% series susceptance
            b = b ./ tap;
            Pfinj = b .* (-branch(:, SHIFT) * pi/180);
            obj.B = sparse( ...
                [1:nl 1:nl nl+1:2*nl nl+1:2*nl]', ...
                [1:nl nl+1:2*nl 1:nl nl+1:2*nl]', ...
                [b; -b; -b; b], ...
                2*nl, 2*nl );
            obj.p = [Pfinj; -Pfinj];
        end

        function add_opf_constraints(obj, nm, om, dm, mpopt)
            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

            %% find branches with flow limits
            mpc = dm.mpc;
            il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);
            nl2 = length(il);         %% number of constrained lines

            if nl2
                %% limits
                flow_max = mpc.branch(il, RATE_A)/mpc.baseMVA;  %% RATE_A

                %% branch flow constraints
                [B, K, p] = obj.get_params(il);
                Af = B * obj.C';
                om.add_lin_constraint('Pf', Af, -p-flow_max, -p+flow_max, ...
                    {nm.va.order(:).name});
            end

            %% branch voltage angle difference limits
            nb = size(mpc.bus, 1);      %% number of buses
            [Aang, lang, uang, iang] = makeAang(mpc.baseMVA, mpc.branch, nb, mpopt);
            om.add_lin_constraint('ang', Aang, lang, uang, {'Va'});
            om.userdata.iang = iang;
        end
    end     %% methods
end         %% classdef
