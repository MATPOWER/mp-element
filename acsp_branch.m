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

        function d2H = apparent_power_hess(obj, x, lam, asm, idx)
            %% branch squared apparent power flow Hessian
            x_ = asm.x2x_(x);           %% convert real to complex x
            nlam = length(lam);

            [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1, idx);
            n = length(S);
            dSc = spdiags(conj(S), 0, n, n);
            dlam = spdiags(lam, 0, nlam, nlam);
            Sx = [Sv1 Sv2 Szr Szi];
            
            d2H = 2 * real( obj.port_inj_power_hess(x_, dSc * lam, 1, idx) + ...
                        Sx.' * dlam * conj(Sx) );
        end

        function d2H = active_power_hess(obj, x, lam, asm, idx)
            %% branch active power flow Hessian
            x_ = asm.x2x_(x);           %% convert real to complex x

            d2H = real( obj.port_inj_power_hess(x_, lam, 1, idx) );
        end

        function d2H = active_power2_hess(obj, x, lam, asm, idx)
            %% branch squared active power flow Hessian
            x_ = asm.x2x_(x);           %% convert real to complex x
            nlam = length(lam);

            [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1, idx);
            n = length(S);
            dP = spdiags(real(S), 0, n, n);
            dlam = spdiags(lam, 0, nlam, nlam);
            Px = real([Sv1 Sv2 Szr Szi]);
            
            d2H = 2 * real( obj.port_inj_power_hess(x_, dP * lam, 1, idx) + ...
                        Px.' * dlam * Px );
        end

        function d2H = current_hess(obj, x, lam, asm, idx)
            %% branch squared current Hessian
            x_ = asm.x2x_(x);           %% convert real to complex x
            nlam = length(lam);

            [I, Iv1, Iv2, Izr, Izi] = obj.port_inj_current(x_, 1, idx);
            n = length(I);
            dIc = spdiags(conj(I), 0, n, n);
            dlam = spdiags(lam, 0, nlam, nlam);
            Ix = [Iv1 Iv2 Izr Izi];
            
            d2H = 2 * real( obj.port_inj_power_hess(x_, dIc * lam, 1, idx) + ...
                        Ix.' * dlam * conj(Ix) );
        end
    end     %% methods
end         %% classdef
