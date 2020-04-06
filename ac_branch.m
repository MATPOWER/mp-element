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

        function [h, dh] = apparent_power_fcn(obj, x, asm, idx, hmax)
            %% branch squared apparent power flow constraints
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
            h = conj(S) .* S - hmax;
        end

        function [h, dh] = active_power_fcn(obj, x, asm, idx, hmax)
            %% branch active power flow constraints
            x_ = asm.x2x_(x);           %% convert real to complex x

            %% get port power injections with derivatives
            if nargout > 1
                [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1, idx);
                P = real(S);
                dh = real([Sv1 Sv2 Szr Szi]);
            else
                P = real(obj.port_inj_power(x, 1, idx));
            end
            h = P - hmax;
        end

        function [h, dh] = active_power2_fcn(obj, x, asm, idx, hmax)
            %% branch squared active power flow constraints
            x_ = asm.x2x_(x);           %% convert real to complex x

            %% get port power injections with derivatives
            if nargout > 1
                [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1, idx);
                Sx = [Sv1 Sv2 Szr Szi];
                P = real(S);
                n = length(S);
                dP = spdiags(P, 0, n, n);
                dh = 2 * real(dP) * real(Sx);
            else
                P = real(obj.port_inj_power(x, 1, idx));
            end
            h = P .* P - hmax;
        end

        function [h, dh] = current_fcn(obj, x, asm, idx, hmax)
            %% branch squared current constraints
            x_ = asm.x2x_(x);           %% convert real to complex x

            %% get port current injections with derivatives
            if nargout > 1
                [I, Iv1, Iv2, Izr, Izi] = obj.port_inj_current(x_, 1, idx);
                Ix = [Iv1 Iv2 Izr Izi];
                n = length(I);
                dI = spdiags(I, 0, n, n);
                dh = 2 * (real(dI) * real(Ix) + imag(dI) * imag(Ix));
            else
                I = obj.port_inj_power(x, 1, idx);
            end
            h = conj(I) .* I - hmax;
        end

        function add_opf_constraints(obj, asm, om, mpc, mpopt)
            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

            %% find branches with flow limits
            il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);
            nl2 = length(il);         %% number of constrained lines

            if nl2
                %% port indexes
                nl = obj.nk;
                idx = [il; nl+il];

                %% limits
                flow_max = mpc.branch(il, RATE_A)/mpc.baseMVA;  %% RATE_A

                %% branch flow constraints
                lim_type = upper(mpopt.opf.flow_lim(1));
                if lim_type == 'S'
                    fcn_flow = @(x)apparent_power_fcn(obj, x, asm, idx, ...
                                                    [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)apparent_power_hess(obj, x, lam, asm, idx);
                elseif lim_type == 'P'
                    fcn_flow = @(x)active_power_fcn(obj, x, asm, idx, ...
                                                    [flow_max; flow_max]);
                    hess_flow = @(x, lam)active_power_hess(obj, x, lam, asm, idx);
                elseif lim_type == '2' || lim_type == 'P'
                    fcn_flow = @(x)active_power2_fcn(obj, x, asm, idx, ...
                                                    [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)active_power2_hess(obj, x, lam, asm, idx);
                elseif lim_type == 'I'
                    fcn_flow = @(x)current_fcn(obj, x, asm, idx, ...
                                                    [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)current_hess(obj, x, lam, asm, idx);
                else
                    error('ac_branch/add_opf_constraints: MPOPT.opf.flow_lim = ''%s'' not yet implemented.', mpopt.opf.flow_lim);
                end
            
                om.add_nln_constraint({'Sf', 'St'}, [nl2;nl2], 0, fcn_flow, hess_flow);
            end

            %% call parent
            add_opf_constraints@mp_branch(obj, asm, om, mpc, mpopt);
        end
    end     %% methods
end         %% classdef
