classdef acc_model < ac_model
%ACC_MODEL  MATPOWER Model class for AC cartesian voltage models.
%   Each concrete MATPOWER Element class must inherit, at least indirectly,
%   from both MP_ELEMENT and MP_MODEL.
%
%   Subclass of AC_MODEL.
%   MP_MODEL provides properties and methods related to the specific
%   model and formulation (e.g. DC version, AC polar power version, etc.)
%
%   Properties
%       (model parameters inherited from AC_MODEL)
%
%   Methods
%       model_name() - returns string w/name of model/formulation ('AC-cartesian model')
%       model_tag() - returns string w/short label for model/formulation ('acc')

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function name = model_name(obj)
            name = 'AC-cartesian model';
        end
        function tag = model_tag(obj)
            tag = 'acc';
        end
        function vtypes = model_vvars(obj)
            vtypes = {'vr', 'vi'};
        end

        function [Iu, Iw] = port_inj_current_jac(obj, ...
                n, v_, Y, M, invdiagvic, diagSlincJ)
            % [Iu, Iw] = obj.port_inj_current_jac(...)

            %% intermediate terms
            E = invdiagvic * (conj(M) - invdiagvic * diagSlincJ);

            %% linear current term
            Ilin_u = Y;
            Ilin_w = 1j * Y;

            %% current from linear power term
            IS_u = E;
            IS_w = -1j * E;

            %% combine
            Iu = Ilin_u + IS_u;
            Iw = Ilin_w + IS_w;
        end

        function [Iuu, Iuw, Iww] = port_inj_current_hess_v(obj, x_, lam, v_, z_, diaginvic, Y, M, diagSlincJ, dlamJ)
            % [Iuu, Iuw, Iww] = obj.port_inj_current_hess_v(x_, lam)
            % [Iuu, Iuw, Iww] = obj.port_inj_current_hess_v(x_, lam, sysx)
            % [Iuu, Iuw, Iww] = obj.port_inj_current_hess_v(x_, lam, sysx, idx)
            % [...] = obj.port_inj_current_hess_vz(x_, lam, v_, z_, diaginvic, Y, M, diagSlincJ, dlamJ)

            if nargin < 10
                sysx = v_;
                if nargin < 5
                    idx = []
                else
                    idx = z_;
                end
                [Y, M, N, s] = obj.get_params(idx, {'Y', 'M', 'N', 's'});
                [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

                %% compute linear power injections
                if isempty(z_)
                    Slin = M*v_ + s;
                else
                    Slin = M*v_ + N*z_ + s;
                end

                n  = length(v_);    %% number of all port voltages
                ni = length(vi_);   %% number of selected port voltages
                diaginvic = sparse(1:ni, 1:ni, 1 ./ conj(vi_), ni, ni);
                if isempty(idx)     %% all ports
                    diagSlincJ = sparse(1:n, 1:n, conj(Slin), n, n);
                    dlamJ      = sparse(1:n, 1:n, lam, n, n);
                else                %% selected ports
                    diagSlincJ = sparse(1:ni, idx, conj(Slin), ni, n);
                    dlamJ      = sparse(1:ni, idx, lam, ni, n);
                end
            else
                n  = length(v_);    %% number of all port voltages
                ni = length(lam);
            end

            %% intermediate terms
            % A = diaginvic * conj(M);
            % B = diaginvic * conj(N);
            % C = diaginvic * diaginvic * diagSlincJ;
            D = (diaginvic * dlamJ).';
            E = diaginvic * (conj(M) - diaginvic * diagSlincJ);
            % E = A - C;
            F = D * E;
            G = -(F.' + F);
            % H = -D * B;

            %% linear current term
            %% second derivatives all zero

            %% current from linear power term
            Iuu = G;
            Iuw = -1j * G;
            Iww = -G;
        end

        function [Iuzr, Iuzi, Iwzr, Iwzi] = port_inj_current_hess_vz(obj, x_, lam, v_, z_, diaginvic, N, dlamJ)
            % [Iuzr, Iuzi, Iwzr, Iwzi] = obj.port_inj_current_hess_vz(x_, lam)
            % [...] = obj.port_inj_current_hess_vz(x_, lam, sysx)
            % [...] = obj.port_inj_current_hess_vz(x_, lam, sysx, idx)
            % [...] = obj.port_inj_current_hess_vz(x_, lam, v_, z_, diaginvic, N, dlamJ)

            if nargin < 8
                sysx = v_;
                if nargin < 5
                    idx = []
                else
                    idx = z_;
                end
                N = obj.get_params(idx, 'N');
                [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

                n  = length(v_);    %% number of all port voltages
                ni = length(vi_);   %% number of selected port voltages
                diaginvic = sparse(1:ni, 1:ni, 1 ./ conj(vi_), ni, ni);
                if isempty(idx)     %% all ports
                    dlamJ = sparse(1:n, 1:n, lam, n, n);
                else                %% selected ports
                    dlamJ = sparse(1:ni, idx, lam, ni, n);
                end
            else
                n  = length(v_);    %% number of all port voltages
                ni = length(lam);
            end

            %% intermediate terms
            % A = diaginvic * conj(M);
            B = diaginvic * conj(N);
            % C = diaginvic * diaginvic * diagSlincJ;
            D = (diaginvic * dlamJ).';
            % E = diaginvic * (conj(M) - diaginvic * diagSlincJ);
            % E = A - C;
            % F = D * E;
            % G = -(F.' + F);
            H = -D * B;

            %% current from linear power term
            Iuzr = H;
            Iuzi = -1j * H;
            Iwzr = Iuzi;
            Iwzi = -H;
        end

        function [Su, Sw] = port_inj_power_jac(obj, ...
                n, v_, Y, M, diagv, diagvi, diagIlincJ)
            % [Su, Sw] = obj.port_inj_power_jac(...)

            %% intermediate terms
            A = diagIlincJ;
            B = diagvi * conj(Y);

            %% linear power term
            Slin_u = M;
            Slin_w = 1j * M;

            %% power from linear current term
            SI_u = A + B;
            SI_w = 1j * (A - B);

            %% combine
            Su = Slin_u + SI_u;
            Sw = Slin_w + SI_w;
        end

        function [Suu, Suw, Sww] = port_inj_power_hess_v(obj, x_, lam, v_, z_, diagvi, Y, M, diagIlincJ, dlamJ)
            % [Suu, Suw, Sww] = obj.port_inj_power_hess_v(x_, lam)
            % [Suu, Suw, Sww] = obj.port_inj_power_hess_v(x_, lam, sysx)
            % [Suu, Suw, Sww] = obj.port_inj_power_hess_v(x_, lam, sysx, idx)
            % [...] = obj.port_inj_power_hess_v(x_, lam, v_, z_, diagvi, Y, M, diagIlincJ, dlamJ)

            if nargin < 10
                sysx = v_;
                if nargin < 5
                    idx = []
                else
                    idx = z_;
                end
                [Y, L, M, i] = obj.get_params(idx, {'Y', 'L', 'M', 'i'});
                [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

                %% compute linear current injections
                if isempty(z_)
                    Ilin = Y*v_ + i;
                else
                    Ilin = Y*v_ + L*z_ + i;
                end

                n  = length(v_);    %% number of all port voltages
                ni = length(vi_);   %% number of selected port voltages
                diagvi = sparse(1:ni, 1:ni, vi_, ni, ni);
                if isempty(idx)     %% all ports
                    diagIlincJ = sparse(1:n, 1:n, conj(Ilin), n, n);
                    dlamJ      = sparse(1:n, 1:n, lam, n, n);
                else                %% selected ports
                    diagIlincJ = sparse(1:ni, idx, conj(Ilin), ni, n);
                    dlamJ      = sparse(1:ni, idx, lam, ni, n);
                end
            else
                n  = length(v_);    %% number of all port voltages
                ni = length(lam);
            end

            %% intermediate terms
            D = dlamJ.';
            E = D * conj(Y);
            F = E + E.';
            G = 1j * (E - E.');

            %% linear power term
            %% second derivatives all zero

            %% power from linear current term
            Suu = F;
            Suw = G.';
            Sww = Suu;
        end

        function [Suzr, Suzi, Swzr, Swzi] = port_inj_power_hess_vz(obj, x_, lam, v_, z_, diagvi, L, dlamJ)
            % [Suzr, Suzi, Swzr, Swzi] = obj.port_inj_power_hess_vz(x_, lam)
            % [...] = obj.port_inj_power_hess_vz(x_, lam, sysx)
            % [...] = obj.port_inj_power_hess_vz(x_, lam, sysx, idx)
            % [...] = obj.port_inj_power_hess_vz(x_, lam, v_, z_, diagvi, L, dlamJ)

            if nargin < 8
                sysx = v_;
                if nargin < 5
                    idx = []
                else
                    idx = z_;
                end
                L = obj.get_params(idx, 'L');
                [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

                n  = length(v_);    %% number of all port voltages
                ni = length(vi_);   %% number of selected port voltages
                diagvi = sparse(1:ni, 1:ni, vi_, ni, ni);
                if isempty(idx)     %% all ports
                    dlamJ = sparse(1:n, 1:n, lam, n, n);
                else                %% selected ports
                    dlamJ = sparse(1:ni, idx, lam, ni, n);
                end
            else
                n  = length(v_);    %% number of all port voltages
                ni = length(lam);
            end

            %% intermediate terms
            D = dlamJ.';
            H = D * conj(L);

            %% power from linear current term
            Suzr = H;
            Suzi = -1j * H;
            Swzr = 1j * H;
            Swzi = H;
        end
    end     %% methods
end         %% classdef
