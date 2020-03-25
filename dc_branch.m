classdef dc_branch < mp_branch & dc_model

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
        %% constructor
        function obj = dc_branch(varargin)
            obj@mp_branch(varargin{:});
        end

        function obj = build_params(obj, asm, mpc)
            build_params@mp_branch(obj, asm, mpc);  %% call parent
            define_constants;
            branch = mpc.branch;
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

        function add_opf_constraints(obj, asm, om, mpc, mpopt)
            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

            %% find branches with flow limits
            il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);
            nl2 = length(il);         %% number of constrained lines

            %% limits
            lims = mpc.branch(il, RATE_A)/mpc.baseMVA;  %% RATE_A

            %% branch flow constraints
            [B, K, p] = obj.get_params(il);
            C = obj.getC();
            Af = B*C';
            om.add_lin_constraint('Pf', Af, -p-lims, -p+lims, {asm.va.order(:).name});

            %% call parent
            add_opf_constraints@mp_branch(obj, asm, om, mpc, mpopt);
        end
    end     %% methods
end         %% classdef
