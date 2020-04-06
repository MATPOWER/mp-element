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
                    fcn_flow = @(x)port_apparent_power_lim_fcn(obj, x, asm, idx, ...
                                                    [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_apparent_power_lim_hess(obj, x, lam, asm, idx);
                elseif lim_type == 'P'
                    fcn_flow = @(x)port_active_power_lim_fcn(obj, x, asm, idx, ...
                                                    [flow_max; flow_max]);
                    hess_flow = @(x, lam)port_active_power_lim_hess(obj, x, lam, asm, idx);
                elseif lim_type == '2' || lim_type == 'P'
                    fcn_flow = @(x)port_active_power2_lim_fcn(obj, x, asm, idx, ...
                                                    [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_active_power2_lim_hess(obj, x, lam, asm, idx);
                elseif lim_type == 'I'
                    fcn_flow = @(x)port_current_lim_fcn(obj, x, asm, idx, ...
                                                    [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_current_lim_hess(obj, x, lam, asm, idx);
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
