classdef ac_model < mp_model
%AC_MODEL  MATPOWER Model base class for AC models.
%   Each concrete MATPOWER Element class must inherit, at least indirectly,
%   from both MP_ELEMENT and MP_MODEL.
%
%   Subclass of MP_MODEL.
%   MP_MODEL provides propoerties and methods related to the specific
%   model and formulation (e.g. DC version, AC polar power version, etc.)
%
%   AC_MODEL defines:
%       linear current injection       = Y v + L z + i
%       linear complex power injection = M v + N z + s
%
%   Properties
%       (model parameters)
%       params - cell array of model parameter field names
%       Y - np*nk x nn matrix
%       L - np*nk x nz matrix
%       M - np*nk x nn matrix
%       N - np*nk x nz matrix
%       i - np*nk x 1 matrix
%       s - np*nk x 1 matrix
%
%   Methods
%       model_params() - cell array of names of model parameters
%                        {'Y', 'L', 'M', 'N', 'i', 's'}

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        %% model parameters
        Y = [];
        L = [];
        M = [];
        N = [];
        i = [];
        s = [];
        param_ncols = struct('Y', 2, 'L', 3, 'M', 2, 'N', 3, 'i', 1, 's', 1);
            %% num of columns for each parameter, where 1 = 1, 2 = np, 3 = nz
    end

    methods
        function name = model_name(obj)
            name = 'AC model';
        end
        function tag = model_tag(obj)
            tag = 'ac';
        end
        function params = model_params(obj)
           params = {'Y', 'L', 'M', 'N', 'i', 's'};
        end
        function vtypes = model_zvars(obj)
            vtypes = {'zr', 'zi'};
        end

        function [I, Iv1, Iv2, Izr, Izi] = port_inj_current(obj, x, sysx, idx)
            % I = obj.port_inj_current(x, sysx)
            % I = obj.port_inj_current(x, sysx, idx)
            % [I, Iv1, Iv2] = obj.port_inj_current(...)
            % [I, Iv1, Iv2, Izr, Izi] = obj.port_inj_current(...)
            % sysx : 1 = system x, 0 = class aggregate x

            if nargin < 4
                idx = [];
                if nargin < 3
                    sysx = 1;
                end
            end
            
            [Y, L, M, N, i, s] = obj.get_params(idx);
            [v, z, vi] = obj.x2vz(x, sysx, idx);

            %% compute linear current injections and power injections
            if isempty(z)
                Slin = M*v + s;
                I = Y*v + i + conj(Slin ./ vi);
            else
                Slin = M*v + N*z + s;
                I = Y*v + L*z + i + conj(Slin ./ vi);
%                 I = Y*v + L*z + i * ones(1, size(x, 2)) + conj(Slin ./ vi);
            end

            if nargout > 1
                %% intermediate terms
                n  = length(v);     %% number of all port voltages
                ni = length(vi);    %% number of selected port voltages
                diagv = sparse(1:n, 1:n, v, n, n);
                invdiagvic = sparse(1:ni, 1:ni, 1./conj(vi), ni, ni);
                if isempty(idx)     %% all ports
                    diagSlincJ = sparse(1:n, 1:n, conj(Slin), n, n);
                else                %% selected ports
                    diagSlincJ = sparse(1:ni, idx, conj(Slin), ni, n);
                end
                
                [Iv1, Iv2] = obj.port_inj_current_jac( ...
                        n, v, Y, M, diagv, invdiagvic, diagSlincJ);

                if nargout >= 4
                    %% linear current term
                    Ilin_zr = L;
                    Ilin_zi = 1j * L;

                    %% current from linear power term
                    IS_zr = invdiagvic * conj(N);   %% C
                    IS_zi = -1j * IS_zr;            %% -jC

                    %% combine
                    Izr = Ilin_zr + IS_zr;
                    Izi = Ilin_zi + IS_zi;
                end

                if sysx
                    Ct = obj.getC('tr');
                    Iv1 = Iv1 * Ct;
                    Iv2 = Iv2 * Ct;
                    if nargout >= 4
                        Dt = obj.getD('tr');
                        Izr = Izr * Dt;
                        Izi = Izi * Dt;
                    end
                end
            end
        end

        function [S, Sv1, Sv2, Szr, Szi] = port_inj_power(obj, x, sysx, idx)
            % S = obj.port_inj_power(x, sysx)
            % S = obj.port_inj_power(x, sysx, idx)
            %   
            % [S, Sv1, Sv2] = obj.port_inj_power(...)
            % [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(...)
            % sysx : 1 = system x, 0 = class aggregate x

            if nargin < 4
                idx = [];
                if nargin < 3
                    sysx = 1;
                end
            end
            
            [Y, L, M, N, i, s] = obj.get_params(idx);
            [v, z, vi] = obj.x2vz(x, sysx, idx);

            %% compute linear current injections and power injections
            if isempty(z)
                Ilin = Y*v + i;
                S = vi .* conj(Ilin) + M*v + s;
            else
                Ilin = Y*v + L*z + i;
                S = vi .* conj(Ilin) + M*v + N*z + s;
%                 S = vi .* conj(Ilin) + M*v + N*z + s * ones(1, size(x, 2));
            end

            if nargout > 1
                %% intermediate terms
                n  = length(v);     %% number of all port voltages
                ni = length(vi);    %% number of selected port voltages
                diagv  = sparse(1:n, 1:n, v, n, n);
                diagvi = sparse(1:ni, 1:ni, vi, ni, ni);
                if isempty(idx)     %% all ports
                    diagIlincJ = sparse(1:n, 1:n, conj(Ilin), n, n);
                else                %% selected ports
                    diagIlincJ = sparse(1:ni, idx, conj(Ilin), ni, n);
                end
                
                [Sv1, Sv2] = obj.port_inj_power_jac( ...
                        n, v, Y, M, diagv, diagvi, diagIlincJ);

                if nargout >= 4
                    %% linear power term
                    Slin_zr = N;
                    Slin_zi = 1j * N;

                    %% power from linear current term
                    SI_zr = diagvi * conj(L);   %% E
                    SI_zi = -1j * SI_zr;        %% -jE

                    %% combine
                    Szr = Slin_zr + SI_zr;
                    Szi = Slin_zi + SI_zi;
                end

                if sysx
                    Ct = obj.getC('tr');
                    Sv1 = Sv1 * Ct;
                    Sv2 = Sv2 * Ct;
                    if nargout >= 4
                        Dt = obj.getD('tr');
                        Szr = Szr * Dt;
                        Szi = Szi * Dt;
                    end
                end
            end
        end

        function [h, dh] = port_apparent_power_lim_fcn(obj, x, asm, idx, hmax)
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

        function [h, dh] = port_active_power_lim_fcn(obj, x, asm, idx, hmax)
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

        function [h, dh] = port_active_power2_lim_fcn(obj, x, asm, idx, hmax)
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

        function [h, dh] = port_current_lim_fcn(obj, x, asm, idx, hmax)
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

        function d2H = port_apparent_power_lim_hess(obj, x, lam, asm, idx)
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

        function d2H = port_active_power_lim_hess(obj, x, lam, asm, idx)
            %% branch active power flow Hessian
            x_ = asm.x2x_(x);           %% convert real to complex x

            d2H = real( obj.port_inj_power_hess(x_, lam, 1, idx) );
        end

        function d2H = port_active_power2_lim_hess(obj, x, lam, asm, idx)
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

        function d2H = port_current_lim_hess(obj, x, lam, asm, idx)
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
