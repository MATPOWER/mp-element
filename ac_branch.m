classdef ac_branch < mp_branch% & ac_model

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
        function obj = ac_branch(varargin)
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

            tap = tap .* exp(1j*pi/180 * branch(:, SHIFT)); %% add phase shifters
            Ys = stat ./ (branch(:, BR_R) + 1j * branch(:, BR_X));  %% series admittance
            Bc = stat .* branch(:, BR_B);   %% line charging susceptance
            Ytt = Ys + 1j*Bc/2;
            Yff = Ytt ./ (tap .* conj(tap));
            Yft = - Ys ./ conj(tap);
            Ytf = - Ys ./ tap;

            obj.Y = sparse( ...
                [1:nl 1:nl nl+1:2*nl nl+1:2*nl]', ...
                [1:nl nl+1:2*nl 1:nl nl+1:2*nl]', ...
                [Yff; Yft; Ytf; Ytt], 2*nl, 2*nl );
        end

        function [h, dh] = opf_branch_flow_fcn(obj, x, asm, idx, lims)
            %% branch power flow constraints
            x_ = asm.x2x_(x);           %% convert real to complex x

            %% get port power injections with derivatives
            if nargout > 1
                [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1, idx);
                Sx = [Sv1 Sv2 Szr Szi];
                n = length(S);
                dS = spdiags(S, 0, n, n);
                dh = 2 * (real(dS) * real(Sx) + imag(dS) * imag(Sx));
            else
                S = obj.port_inj_power(x, 1, idx);
            end
            h = conj(S) .* S - lims;
        end

        function add_opf_constraints(obj, asm, om, mpc, mpopt)
        
            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
            
            %% find branches with flow limits
            il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);
            nl2 = length(il);         %% number of constrained lines

            %% port indexes
            nl = obj.nk;
            idx = [il; nl+il];

            %% limits
            lims = (mpc.branch(il, RATE_A)/mpc.baseMVA) .^ 2;   %% square of RATE_A

            %% branch flow constraints
            fcn_flow = @(x)opf_branch_flow_fcn(obj, x, asm, idx, [lims; lims]);
            hess_flow = @(x, lam)opf_branch_flow_hess(obj, x, lam, asm, idx);
            om.add_nln_constraint({'Sf', 'St'}, [nl2;nl2], 0, fcn_flow, hess_flow);
        end
    end     %% methods
end         %% classdef
