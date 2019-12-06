classdef ac_branch < mp_branch & acsp_model

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
    end     %% methods
end         %% classdef
