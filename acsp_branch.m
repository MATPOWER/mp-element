classdef acsp_branch < ac_branch & acsp_model

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
        function obj = acsp_branch(varargin)
            obj@ac_branch(varargin{:});
        end

        function d2H = opf_branch_flow_hess(obj, x, lam, asm, idx)
            %% branch power flow Hessian
            x_ = asm.x2x_(x);           %% convert real to complex x
            nlam = length(lam);

            [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1, idx);
            n = length(S);
            dSc = spdiags(conj(S), 0, n, n);
            dlam = spdiags(lam, 0, nlam, nlam);
            Sx = [Sv1 Sv2 Szr Szi];
            
            d2H = 2 * real( obj.port_inj_power_hess(x_, dSc * lam, 1, idx) + ...
                        Sx.' * dlam * conj(Sx)  );
        end
    end     %% methods
end         %% classdef
