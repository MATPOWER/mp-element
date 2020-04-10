classdef acp_model < ac_model
%ACP_MODEL  MATPOWER Model class for AC polar voltage models.
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
%       model_name() - returns string w/name of model/formulation ('AC-polar model')
%       model_tag() - returns string w/short label for model/formulation ('acp')

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
            name = 'AC-polar model';
        end
        function tag = model_tag(obj)
            tag = 'acp';
        end
        function vtypes = model_vvars(obj)
            vtypes = {'va', 'vm'};
        end

        function [Iva, Ivm] = port_inj_current_jac(obj, ...
                n, v_, Y, M, invdiagvic, diagSlincJ)
            % [Iva, Ivm] = obj.port_inj_current_jac(...)

            %% intermediate terms
            diagv = sparse(1:n, 1:n, v_, n, n);
            C = invdiagvic * (diagSlincJ - conj(M * diagv));
            D = sparse(1:n, 1:n, 1 ./ abs(v_), n, n);

            %% linear current term
            Ilin_vm = Y * diagv;
            Ilin_va = 1j * Ilin_vm;
            Ilin_vm = Ilin_vm * D;

            %% current from linear power term
            IS_va = 1j * C;
            IS_vm = -C * D;

            %% combine
            Iva = Ilin_va + IS_va;
            Ivm = Ilin_vm + IS_vm;
        end

        function H = port_inj_current_hess(obj, x_, lam, sysx, idx)
            % H = obj.port_inj_current_hess(x_, lam)
            % H = obj.port_inj_current_hess(x_, lam, sysx)
            % H = obj.port_inj_current_hess(x_, lam, sysx, idx)

            if nargin < 5
                idx = [];
                if nargin < 4
                    sysx = 1;
                end
            end

            [Y, M, N, s] = obj.get_params(idx, {'Y', 'M', 'N', 's'});
            [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

            %% compute linear power injections
            if isempty(z_)
                Slin = M*v_ + s;
            else
                Slin = M*v_ + N*z_ + s;
            end

            %% intermediate terms
            nz = length(z_);
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

            [Ivava, Ivavm, Ivmvm] = ...
                obj.port_inj_current_hess_v(x_, lam, v_, z_, diaginvic, Y, M, diagSlincJ, dlamJ);
            [Ivazr, Ivazi, Ivmzr, Ivmzi] = ...
                obj.port_inj_current_hess_vz(x_, lam, v_, z_, diaginvic, N, dlamJ);

            H = [   Ivava   Ivavm   Ivazr   Ivazi;
                    Ivavm.' Ivmvm   Ivmzr   Ivmzi;
                    Ivazr.' Ivmzr.' sparse(nz, 2*nz);
                    Ivazi.' Ivmzi.' sparse(nz, 2*nz)  ];

            %% convert for system x_, if necessary
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

            %% general nonlinear current
            if ~isempty(obj.inln_hess)
                H = H + obj.inln_hess(x_, lam, sysx, idx);
            elseif ~isempty(obj.snln_hess)
                error('Nonlinear current Hessian not defined for corresponding nonlinear power Hessian.')
            end
        end

        function [Ivava, Ivavm, Ivmvm] = port_inj_current_hess_v(obj, x_, lam, v_, z_, diaginvic, Y, M, diagSlincJ, dlamJ)
            % [Ivava, Ivavm, Ivmvm] = obj.port_inj_current_hess_v(x_, lam)
            % [Ivava, Ivavm, Ivmvm] = obj.port_inj_current_hess_v(x_, lam, sysx)
            % [Ivava, Ivavm, Ivmvm] = obj.port_inj_current_hess_v(x_, lam, sysx, idx)
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
            diagv  = sparse(1:n, 1:n, v_, n, n);
            A = diaginvic * diagSlincJ;
            B = diaginvic * conj(M);
            % B2 = B * conj(diagv);
            % C = A - B2;
            D = sparse(1:n, 1:n, 1 ./ abs(v_), n, n);
            % E = diaginvic * conj(N);
            F = dlamJ.';
            G = F * B * conj(diagv);    %% F * B2
            dBtlam = sparse(1:n, 1:n, B.' * lam, n, n);
            H = dBtlam * conj(diagv);
            K = (F * A).';
            GG = G + G.';
            LL = GG - H - K;
            MM = D * (2*K - GG) * D;

            %% linear current term
            dYtlam = sparse(1:n, 1:n, Y.' * lam, n, n);
            Ilin_vava = -dYtlam * diagv;
            Ilin_vavm = -j * Ilin_vava * D;

            %% current from linear power term
            IS_vava = LL;
            IS_vavm = 1j * LL.' * D;
            IS_vmvm = MM;

            %% combine
            Ivava = Ilin_vava + IS_vava;
            Ivavm = Ilin_vavm + IS_vavm;
            Ivmvm = IS_vmvm;
        end

        function [Ivazr, Ivazi, Ivmzr, Ivmzi] = port_inj_current_hess_vz(obj, x_, lam, v_, z_, diaginvic, N, dlamJ)
            % [Ivazr, Ivazi, Ivmzr, Ivmzi] = obj.port_inj_current_hess_vz(x_, lam)
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
            D = sparse(1:n, 1:n, 1 ./ abs(v_), n, n);
            E = diaginvic * conj(N);
            NN = dlamJ.' * E;

            %% current from linear power term
            Ivazr = 1j * NN;
            Ivazi = NN;
            Ivmzr = -D * NN;
            Ivmzi = -1j * Ivmzr;
        end

        function [Sva, Svm] = port_inj_power_jac(obj, ...
                n, v_, Y, M, diagv, diagvi, diagIlincJ)
            % [Sva, Svm] = obj.port_inj_power_jac(...)

            %% intermediate terms
            A = diagvi * diagIlincJ;
            B = diagvi * conj(Y);
            C = B * conj(diagv);
            D = sparse(1:n, 1:n, 1 ./ abs(v_), n, n);

            %% linear power term
            Slin_vm = M * diagv;
            Slin_va = 1j * Slin_vm;
            Slin_vm = Slin_vm * D;

            %% power from linear current term
            SI_va = 1j * (A - C);
            SI_vm = (A + C) * D;

            %% combine
            Sva = Slin_va + SI_va;
            Svm = Slin_vm + SI_vm;
        end

        function H = port_inj_power_hess(obj, x_, lam, sysx, idx)
            % H = obj.port_inj_power_hess(x_, lam)
            % H = obj.port_inj_power_hess(x_, lam, sysx)
            % H = obj.port_inj_power_hess(x_, lam, sysx, idx)

            if nargin < 5
                idx = [];
                if nargin < 4
                    sysx = 1;
                end
            end

            [Y, L, M, i] = obj.get_params(idx, {'Y', 'L', 'M', 'i'});
            [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);

            %% compute linear current injections
            if isempty(z_)
                Ilin = Y*v_ + i;
            else
                Ilin = Y*v_ + L*z_ + i;
            end

            %% intermediate terms
            nz = length(z_);
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

            [Svava, Svavm, Svmvm] = ...
                obj.port_inj_power_hess_v(x_, lam, v_, z_, diagvi, Y, M, diagIlincJ, dlamJ);
            [Svazr, Svazi, Svmzr, Svmzi] = ...
                obj.port_inj_power_hess_vz(x_, lam, v_, z_, diagvi, L, dlamJ);

            H = [   Svava   Svavm   Svazr   Svazi;
                    Svavm.' Svmvm   Svmzr   Svmzi;
                    Svazr.' Svmzr.' sparse(nz, 2*nz);
                    Svazi.' Svmzi.' sparse(nz, 2*nz)  ];

            %% convert for system x_, if necessary
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

            %% general nonlinear power
            if ~isempty(obj.snln_hess)
                H = H + obj.snln_hess(x_, lam, sysx, idx);
            elseif ~isempty(obj.inln_hess)
                error('Nonlinear power Hessian not defined for corresponding nonlinear current Hessian.')
            end
        end

        function [Svava, Svavm, Svmvm] = port_inj_power_hess_v(obj, x_, lam, v_, z_, diagvi, Y, M, diagIlincJ, dlamJ)
            % [Svava, Svavm, Svmvm] = obj.port_inj_power_hess_v(x_, lam)
            % [Svava, Svavm, Svmvm] = obj.port_inj_power_hess_v(x_, lam, sysx)
            % [Svava, Svavm, Svmvm] = obj.port_inj_power_hess_v(x_, lam, sysx, idx)
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
            diagv  = sparse(1:n, 1:n, v_, n, n);

            A = diagvi * diagIlincJ;
            B = diagvi * conj(Y);
            C = B * conj(diagv);
            D = sparse(1:n, 1:n, 1 ./ abs(v_), n, n);
            % E = diagvi * conj(L);
            F = dlamJ.';
            G = F * C;
            dBtlam = sparse(1:n, 1:n, B.' * lam, n, n);
            H = conj(diagv) * ((F*B).' - dBtlam);
            K = F * (C - A);

            %% linear power term
            dMtlam = sparse(1:n, 1:n, M.' * lam, n, n);
            Slin_vava = -dMtlam * diagv;
            Slin_vavm = -j * Slin_vava * D;

            %% power from linear current term
            SI_vava = H + K;
            SI_vavm = 1j * (H - K).' * D;
            SI_vmvm =  D * (G + G.') * D;

            %% combine
            Svava = Slin_vava + SI_vava;
            Svavm = Slin_vavm + SI_vavm;
            Svmvm = SI_vmvm;
        end

        function [Svazr, Svazi, Svmzr, Svmzi] = port_inj_power_hess_vz(obj, x_, lam, v_, z_, diagvi, L, dlamJ)
            % [Svazr, Svazi, Svmzr, Svmzi] = obj.port_inj_power_hess_vz(x_, lam)
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
            D = sparse(1:n, 1:n, 1 ./ abs(v_), n, n);
            E = diagvi * conj(L);
            LL = dlamJ.' * E;
            M = D * LL;

            %% power from linear current term
            Svazr = 1j * LL;
            Svazi = LL;
            Svmzr = M;
            Svmzi = -1j * M;
        end
    end     %% methods
end         %% classdef
