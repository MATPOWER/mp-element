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
    end     %% methods
end         %% classdef
