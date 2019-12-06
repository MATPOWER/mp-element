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
    end

    methods
%         function name = model_name(obj)
%             name = 'AC model';
%         end
%         function tag = model_tag(obj)
%             tag = 'ac';
%         end
        function params = model_params(obj)
           params = {'Y', 'L', 'M', 'N', 'i', 's'};
        end
        function vtypes = model_zvars(obj)
            vtypes = {'zr', 'zi'};
        end

        function [Y, L, M, N, i, s] = get_params(obj, idx)
            np = obj.nk * obj.np;
            nz = obj.nk * obj.nz;
            if nargin >= 2 && ~isempty(idx) %% selected ports
                ni = length(idx);
                if isempty(obj.Y)
                    Y = sparse(ni, np);
                else
                    Y = obj.Y(idx, :);
                end
                if nargout > 1
                    if isempty(obj.L)
                        L = sparse(ni, nz);
                    else
                        L = obj.L(idx, :);
                    end
                    if nargout > 2
                        if isempty(obj.M)
                            M = sparse(ni, np);
                        else
                            M = obj.M(idx, :);
                        end
                        if nargout > 3
                            if isempty(obj.N)
                                N = sparse(ni, nz);
                            else
                                N = obj.N(idx, :);
                            end
                            if nargout > 4
                                if isempty(obj.i)
                                    i = zeros(ni, 1);
                                else
                                    i = obj.i(idx);
                                end
                                if nargout > 5
                                    if isempty(obj.s)
                                        s = zeros(ni, 1);
                                    else
                                        s = obj.s(idx);
                                    end
                                end
                            end
                        end
                    end
                end
            else                            %% all ports
                if isempty(obj.Y)
                    Y = sparse(np, np);
                else
                    Y = obj.Y;
                end
                if nargout > 1
                    if isempty(obj.L)
                        L = sparse(np, nz);
                    else
                        L = obj.L;
                    end
                    if nargout > 2
                        if isempty(obj.M)
                            M = sparse(np, np);
                        else
                            M = obj.M;
                        end
                        if nargout > 3
                            if isempty(obj.N)
                                N = sparse(np, nz);
                            else
                                N = obj.N;
                            end
                            if nargout > 4
                                if isempty(obj.i)
                                    i = zeros(np, 1);
                                else
                                    i = obj.i;
                                end
                                if nargout > 5
                                    if isempty(obj.s)
                                        s = zeros(np, 1);
                                    else
                                        s = obj.s;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

%         function I = port_inj_current_old(obj, x, sysx, idx)
%             
%             [Y, L, M, N, i, s] = obj.get_params();
%             [v, z] = obj.x2vz(x, sysx);
%             
%             if sysx
%                 Ct = obj.getC('tr');
%                 Dt = obj.getD('tr');
%                 if nargin < 4       %% all ports
%                     I = Y*Ct*v + L*Dt*z + i + conj((M*Ct*v + N*Dt*z + s) ./ (Ct*v));
%                 else                %% selected ports
%                     I = Y(idx, :)*Ct*v + L(idx, :)*Dt*z + i(idx) + ...
%                         conj((M(idx, :)*Ct*v + N(idx, :)*Dt*z + s(idx)) ./ (Ct(idx, :)*v));
%                 end
%             else
%                 if nargin < 4       %% all ports
%                     I = Y*v + L*z + i + conj((M*v + N*z + s) ./ v);
%                 else                %% selected ports
%                     I = Y(idx, :)*v + L(idx, :)*z + i(idx) + ...
%                         conj((M(idx, :)*v + N(idx, :)*z + s(idx)) ./ v(idx));
%                 end
%             end
%         end

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
            Slin = M*v + N*z + s;
            I = Y*v + L*z + i + conj(Slin ./ vi);

            if nargout > 1
                %% intermediate terms
                n  = length(v);     %% number of all port voltages
                ni = length(vi);    %% number of selected port voltages
                diagv = sparse(1:n, 1:n, v, n, n);
                invdiagvic = sparse(1:ni, 1:ni, 1./conj(vi), ni, ni);
                if isempty(idx)     %% all ports
                    diagSlinc = sparse(1:n, 1:n, conj(Slin), n, n);
                else                %% selected ports
                    diagSlinc = sparse(1:ni, idx, conj(Slin), ni, n);
                end
                
                [Iv1, Iv2] = obj.port_inj_current_jac( ...
                        n, v, Y, M, diagv, invdiagvic, diagSlinc);

                if nargout >= 4
                    %% linear current term
                    Ilin_zr = L;
                    Ilin_zi = 1j * L;

                    %% current from linear power term
                    IS_zr = invdiagvic * conj(N);   %% B
                    IS_zi = 1j * IS_zr;             %% jB

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

%         function S = port_inj_power_old(obj, x, sysx, idx)
%             
%             [Y, L, M, N, i, s] = obj.get_params();
%             [v, z] = obj.x2vz(x, sysx);
%             
%             if sysx
%                 Ct = obj.getC('tr');
%                 Dt = obj.getD('tr');
%                 if nargin < 4       %% all ports
%                     S = Ct*v .* conj(Y*Ct*v + L*Dt*z + i) + M*Ct*v + N*Dt*z + s;
%                 else                %% selected ports
%                     S = Ct(idx, :)*v .* conj(Y(idx, :)*Ct*v + L(idx, :)*Dt*z + i(idx)) + ...
%                         M(idx, :)*Ct*v + N(idx, :)*Dt*z + s(idx);
%                 end
%                 if nargout > 1
%                 end
%             else
%                 if nargin < 4       %% all ports
%                     S = v .* conj(Y*v + L*z + i) + M*v + N*z + s;
%                 else                %% selected ports
%                     S = v(idx) .* conj(Y(idx, :)*v + L(idx, :)*z + i) + ...
%                         M(idx, :)*v + N(idx, :)*z + s(idx);
%                 end
%                 if nargout > 1
%                 end
%             end
%         end

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
            Ilin = Y*v + L*z + i;
            S = vi .* conj(Ilin) + M*v + N*z + s;

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
                    SI_zi = 1j * SI_zr;         %% jE

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

        function H = port_inj_power_hess(obj, x, lam, sysx, idx)
            % H = obj.port_inj_power_hess(x, lam)
            % H = obj.port_inj_power_hess(x, lam, sysx)
            % H = obj.port_inj_power_hess(x, lam, sysx, idx)

            if nargin < 5
                idx = [];
            end

            [Y, L, M, N, i] = obj.get_params(idx);
            [v, z, vi] = obj.x2vz(x, sysx, idx);

            %% compute linear current injections
            Ilin = Y*v + L*z + i;

            %% intermediate terms
            nz = length(z);
            n  = length(v);     %% number of all port voltages
            ni = length(vi);    %% number of selected port voltages
            diagvi = sparse(1:ni, 1:ni, vi, ni, ni);
            if isempty(idx)     %% all ports
                diagIlincJ = sparse(1:n, 1:n, conj(Ilin), n, n);
                dlamJ      = sparse(1:n, 1:n, lam, n, n);
            else                %% selected ports
                diagIlincJ = sparse(1:ni, idx, conj(Ilin), ni, n);
                dlamJ      = sparse(1:ni, idx, lam, ni, n);
            end

            [Svava, Svavm, Svmvm] = ...
                obj.port_inj_power_hess_v(x, lam, v, z, diagvi, Y, L, M, diagIlincJ, dlamJ);
            [Svazr, Svazi, Svmzr, Svmzi] = ...
                obj.port_inj_power_hess_vz(x, lam, v, z, diagvi, L, dlamJ);

            H = [   Svava   Svavm   Svazr   Svazi;
                    Svavm.' Svmvm   Svmzr   Svmzi;
                    Svazr.' Svmzr.' sparse(nz, 2*nz);
                    Svazi.' Svmzi.' sparse(nz, 2*nz)  ];
        end
    end     %% methods
end         %% classdef
