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
                n, v, Y, M, invdiagvic, diagSlincJ)
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

        function H = port_inj_current_hess(obj, x, lam, sysx, idx)
            % H = obj.port_inj_current_hess(x, lam)
            % H = obj.port_inj_current_hess(x, lam, sysx)
            % H = obj.port_inj_current_hess(x, lam, sysx, idx)

            if nargin < 5
                idx = [];
                if nargin < 4
                    sysx = 1;
                end
            end

            [Y, M, N, s] = obj.get_params(idx, {'Y', 'M', 'N', 's'});
            [v, z, vi] = obj.x2vz(x, sysx, idx);

            %% compute linear power injections
            if isempty(z)
                Slin = M*v + s;
            else
                Slin = M*v + N*z + s;
            end

            %% intermediate terms
            nz = length(z);
            n  = length(v);     %% number of all port voltages
            ni = length(vi);    %% number of selected port voltages
            diaginvic = sparse(1:ni, 1:ni, 1 ./ conj(vi), ni, ni);
            if isempty(idx)     %% all ports
                diagSlincJ = sparse(1:n, 1:n, conj(Slin), n, n);
                dlamJ      = sparse(1:n, 1:n, lam, n, n);
            else                %% selected ports
                diagSlincJ = sparse(1:ni, idx, conj(Slin), ni, n);
                dlamJ      = sparse(1:ni, idx, lam, ni, n);
            end

            [Iuu, Iuw, Iww] = ...
                obj.port_inj_current_hess_v(x, lam, v, z, diaginvic, Y, M, diagSlincJ, dlamJ);
            [Iuzr, Iuzi, Iwzr, Iwzi] = ...
                obj.port_inj_current_hess_vz(x, lam, v, z, diaginvic, N, dlamJ);

            H = [   Iuu   Iuw   Iuzr   Iuzi;
                    Iuw.' Iww   Iwzr   Iwzi;
                    Iuzr.' Iwzr.' sparse(nz, 2*nz);
                    Iuzi.' Iwzi.' sparse(nz, 2*nz)  ];

            %% convert for system x, if necessary
            if sysx
                C = obj.getC();
                D = obj.getD();
                [mc, nc] = size(C);
                [md, nd] = size(D);
                Ap = [  C sparse(mc, nc+2*nd);
                        sparse(mc,nc) C sparse(mc, 2*nd);
                        sparse(md, 2*nc) D sparse(md,nd);
                        sparse(md, 2*nc+nd) D ];
                H = Ap * H * Ap.';
            end
        end

        function [Iuu, Iuw, Iww] = port_inj_current_hess_v(obj, x, lam, v, z, diaginvic, Y, M, diagSlincJ, dlamJ)
            % [Iuu, Iuw, Iww] = obj.port_inj_current_hess_v(x, lam)
            % [Iuu, Iuw, Iww] = obj.port_inj_current_hess_v(x, lam, sysx)
            % [Iuu, Iuw, Iww] = obj.port_inj_current_hess_v(x, lam, sysx, idx)
            % [...] = obj.port_inj_current_hess_vz(x, lam, v, z, diaginvic, Y, M, diagSlincJ, dlamJ)

            if nargin < 10
                sysx = v;
                if nargin < 5
                    idx = []
                else
                    idx = z;
                end
                [Y, M, N, s] = obj.get_params(idx, {'Y', 'M', 'N', 's'});
                [v, z, vi] = obj.x2vz(x, sysx, idx);

                %% compute linear power injections
                if isempty(z)
                    Slin = M*v + s;
                else
                    Slin = M*v + N*z + s;
                end

                n  = length(v);     %% number of all port voltages
                ni = length(vi);    %% number of selected port voltages
                diaginvic = sparse(1:ni, 1:ni, 1 ./ conj(vi), ni, ni);
                if isempty(idx)     %% all ports
                    diagSlincJ = sparse(1:n, 1:n, conj(Slin), n, n);
                    dlamJ      = sparse(1:n, 1:n, lam, n, n);
                else                %% selected ports
                    diagSlincJ = sparse(1:ni, idx, conj(Slin), ni, n);
                    dlamJ      = sparse(1:ni, idx, lam, ni, n);
                end
            else
                n  = length(v);     %% number of all port voltages
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

        function [Iuzr, Iuzi, Iwzr, Iwzi] = port_inj_current_hess_vz(obj, x, lam, v, z, diaginvic, N, dlamJ)
            % [Iuzr, Iuzi, Iwzr, Iwzi] = obj.port_inj_current_hess_vz(x, lam)
            % [...] = obj.port_inj_current_hess_vz(x, lam, sysx)
            % [...] = obj.port_inj_current_hess_vz(x, lam, sysx, idx)
            % [...] = obj.port_inj_current_hess_vz(x, lam, v, z, diaginvic, N, dlamJ)

            if nargin < 8
                sysx = v;
                if nargin < 5
                    idx = []
                else
                    idx = z;
                end
                N = obj.get_params(idx, 'N');
                [v, z, vi] = obj.x2vz(x, sysx, idx);

                n  = length(v);     %% number of all port voltages
                ni = length(vi);    %% number of selected port voltages
                diaginvic = sparse(1:ni, 1:ni, 1 ./ conj(vi), ni, ni);
                if isempty(idx)     %% all ports
                    dlamJ = sparse(1:n, 1:n, lam, n, n);
                else                %% selected ports
                    dlamJ = sparse(1:ni, idx, lam, ni, n);
                end
            else
                n  = length(v);     %% number of all port voltages
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
                n, v, Y, M, diagv, diagvi, diagIlincJ)
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

        function H = port_inj_power_hess(obj, x, lam, sysx, idx)
            % H = obj.port_inj_power_hess(x, lam)
            % H = obj.port_inj_power_hess(x, lam, sysx)
            % H = obj.port_inj_power_hess(x, lam, sysx, idx)

            if nargin < 5
                idx = [];
                if nargin < 4
                    sysx = 1;
                end
            end

            [Y, L, M, i] = obj.get_params(idx, {'Y', 'L', 'M', 'i'});
            [v, z, vi] = obj.x2vz(x, sysx, idx);

            %% compute linear current injections
            if isempty(z)
                Ilin = Y*v + i;
            else
                Ilin = Y*v + L*z + i;
            end

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

            [Suu, Suw, Sww] = ...
                obj.port_inj_power_hess_v(x, lam, v, z, diagvi, Y, M, diagIlincJ, dlamJ);
            [Suzr, Suzi, Swzr, Swzi] = ...
                obj.port_inj_power_hess_vz(x, lam, v, z, diagvi, L, dlamJ);

            H = [   Suu   Suw   Suzr   Suzi;
                    Suw.' Sww   Swzr   Swzi;
                    Suzr.' Swzr.' sparse(nz, 2*nz);
                    Suzi.' Swzi.' sparse(nz, 2*nz)  ];

            %% convert for system x, if necessary
            if sysx
                C = obj.getC();
                D = obj.getD();
                [mc, nc] = size(C);
                [md, nd] = size(D);
                Ap = [  C sparse(mc, nc+2*nd);
                        sparse(mc,nc) C sparse(mc, 2*nd);
                        sparse(md, 2*nc) D sparse(md,nd);
                        sparse(md, 2*nc+nd) D ];
                H = Ap * H * Ap.';
            end
        end

        function [Suu, Suw, Sww] = port_inj_power_hess_v(obj, x, lam, v, z, diagvi, Y, M, diagIlincJ, dlamJ)
            % [Suu, Suw, Sww] = obj.port_inj_power_hess_v(x, lam)
            % [Suu, Suw, Sww] = obj.port_inj_power_hess_v(x, lam, sysx)
            % [Suu, Suw, Sww] = obj.port_inj_power_hess_v(x, lam, sysx, idx)
            % [...] = obj.port_inj_power_hess_v(x, lam, v, z, diagvi, Y, M, diagIlincJ, dlamJ)

            if nargin < 10
                sysx = v;
                if nargin < 5
                    idx = []
                else
                    idx = z;
                end
                [Y, L, M, i] = obj.get_params(idx, {'Y', 'L', 'M', 'i'});
                [v, z, vi] = obj.x2vz(x, sysx, idx);

                %% compute linear current injections
                if isempty(z)
                    Ilin = Y*v + i;
                else
                    Ilin = Y*v + L*z + i;
                end

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
            else
                n  = length(v);     %% number of all port voltages
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

        function [Suzr, Suzi, Swzr, Swzi] = port_inj_power_hess_vz(obj, x, lam, v, z, diagvi, L, dlamJ)
            % [Suzr, Suzi, Swzr, Swzi] = obj.port_inj_power_hess_vz(x, lam)
            % [...] = obj.port_inj_power_hess_vz(x, lam, sysx)
            % [...] = obj.port_inj_power_hess_vz(x, lam, sysx, idx)
            % [...] = obj.port_inj_power_hess_vz(x, lam, v, z, diagvi, L, dlamJ)

            if nargin < 8
                sysx = v;
                if nargin < 5
                    idx = []
                else
                    idx = z;
                end
                L = obj.get_params(idx, 'L');
                [v, z, vi] = obj.x2vz(x, sysx, idx);

                n  = length(v);     %% number of all port voltages
                ni = length(vi);    %% number of selected port voltages
                diagvi = sparse(1:ni, 1:ni, vi, ni, ni);
                if isempty(idx)     %% all ports
                    dlamJ = sparse(1:n, 1:n, lam, n, n);
                else                %% selected ports
                    dlamJ = sparse(1:ni, idx, lam, ni, n);
                end
            else
                n  = length(v);     %% number of all port voltages
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
