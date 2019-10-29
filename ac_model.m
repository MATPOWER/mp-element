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

        function [Y, L, M, N, i, s] = get_params(obj)
            np = obj.nk * obj.np;
            nz = obj.nk * obj.nz;
            if isempty(obj.Y)
                Y = sparse(np, np);
            else
                Y = obj.Y;
            end
            if isempty(obj.L)
                L = sparse(np, nz);
            else
                L = obj.L;
            end
            if isempty(obj.M)
                M = sparse(np, np);
            else
                M = obj.M;
            end
            if isempty(obj.N)
                N = sparse(np, nz);
            else
                N = obj.N;
            end
            if isempty(obj.i)
                i = zeros(np, 1);
            else
                i = obj.i;
            end
            if isempty(obj.s)
                s = zeros(np, 1);
            else
                s = obj.s;
            end
        end

        function I = port_inj_current(obj, x, sysx, idx)
            % sys x : 1 = system x, 0 = class aggregate x
            
            [Y, L, M, N, i, s] = obj.get_params();
            [v, z] = obj.x2vz(x, sysx);
            
            if sysx
                Ct = horzcat(obj.C{:}).';
                Dt = horzcat(obj.D{:}).';
                if nargin < 4       %% all ports
                    I = Y*Ct*v + L*Dt*z + i + conj((M*Ct*v + N*Dt*z + s) ./ (Ct*v));
                else                %% selected ports
                    I = Y(idx, :)*Ct*v + L(idx, :)*Dt*z + i(idx) + ...
                        conj((M(idx, :)*Ct*v + N(idx, :)*Dt*z + s(idx)) ./ (Ct(idx, :)*v));
                end
            else
                if nargin < 4       %% all ports
                    I = Y*v + L*z + i + conj((M*v + N*z + s) ./ v);
                else                %% selected ports
                    I = Y(idx, :)*v + L(idx, :)*z + i(idx) + ...
                        conj((M(idx, :)*v + N(idx, :)*z + s(idx)) ./ v(idx));
                end
            end
        end

        function S = port_inj_power(obj, x, sysx, idx)
            % sys x : 1 = system x, 0 = class aggregate x
            
            [Y, L, M, N, i, s] = obj.get_params();
            [v, z] = obj.x2vz(x, sysx);
            
            if sysx
                Ct = horzcat(obj.C{:}).';
                Dt = horzcat(obj.D{:}).';
                if nargin < 4       %% all ports
                    S = Ct*v .* conj(Y*Ct*v + L*Dt*z + i) + M*Ct*v + N*Dt*z + s;
                else                %% selected ports
                    S = Ct(idx, :)*v .* conj(Y(idx, :)*Ct*v + L(idx, :)*Dt*z + i(idx)) + ...
                        M(idx, :)*Ct*v + N(idx, :)*Dt*z + s(idx);
                end
            else
                if nargin < 4       %% all ports
                    S = v .* conj(Y*v + L*z + i) + M*v + N*z + s;
                else                %% selected ports
                    S = v(idx) .* conj(Y(idx, :)*v + L(idx, :)*z + i) + ...
                        M(idx, :)*v + N(idx, :)*z + s(idx);
                end
            end
        end
    end     %% methods
end         %% classdef
