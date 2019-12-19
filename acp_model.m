classdef acp_model < ac_model
%ACP_MODEL  MATPOWER Model class for AC polar models.
%   Each concrete MATPOWER Element class must inherit, at least indirectly,
%   from both MP_ELEMENT and MP_MODEL.
%
%   Subclass of AC_MODEL.
%   MP_MODEL provides propoerties and methods related to the specific
%   model and formulation (e.g. DC version, AC polar power version, etc.)
%
%   Properties
%       (model parameters inherited from AC_MODEL)
%
%   Methods
%       model_name() - returns string w/name of model/formulation ('AC-power-polar model')
%       model_tag() - returns string w/short label for model/formulation ('acsp')

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function vtypes = model_vvars(obj)
            vtypes = {'va', 'vm'};
        end

        function [Iva, Ivm] = port_inj_current_jac(obj, ...
                n, v, Y, M, diagv, invdiagvic, diagSlincJ)
            % [Iva, Ivm] = obj.port_inj_current_jac(...)

            %% intermediate terms
            C = invdiagvic * (diagSlincJ - conj(M * diagv));
            D = sparse(1:n, 1:n, 1 ./ abs(v), n, n);

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

        function [Sva, Svm] = port_inj_power_jac(obj, ...
                n, v, Y, M, diagv, diagvi, diagIlincJ)
            % [Sva, Svm] = obj.port_inj_power_jac(...)

            %% intermediate terms
            A = diagvi * diagIlincJ;
            B = diagvi * conj(Y);
            C = B * conj(diagv);
            D = sparse(1:n, 1:n, 1 ./ abs(v), n, n);

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
                J          = sparse(1:n, 1:n, 1, n, n);
            else                %% selected ports
                diagIlincJ = sparse(1:ni, idx, conj(Ilin), ni, n);
                dlamJ      = sparse(1:ni, idx, lam, ni, n);
                J          = sparse(1:ni, idx, 1, ni, n);
            end

            [Svava, Svavm, Svmvm] = ...
                obj.port_inj_power_hess_v(x, lam, v, z, diagvi, Y, L, M, J, diagIlincJ, dlamJ);
            [Svazr, Svazi, Svmzr, Svmzi] = ...
                obj.port_inj_power_hess_vz(x, lam, v, z, diagvi, L, dlamJ);

            H = [   Svava   Svavm   Svazr   Svazi;
                    Svavm.' Svmvm   Svmzr   Svmzi;
                    Svazr.' Svmzr.' sparse(nz, 2*nz);
                    Svazi.' Svmzi.' sparse(nz, 2*nz)  ];

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

        function [Svava, Svavm, Svmvm] = port_inj_power_hess_v(obj, x, lam, v, z, diagvi, Y, L, M, J, diagIlincJ, dlamJ)
            % [Svava, Svavm, Svmvm] = obj.port_inj_power_hess_v(x, lam)
            % [Svava, Svavm, Svmvm] = obj.port_inj_power_hess_v(x, lam, sysx)
            % [Svava, Svavm, Svmvm] = obj.port_inj_power_hess_v(x, lam, sysx, idx)
            % [...] = obj.port_inj_power_hess_vz(x, lam, v, z, diagvi, Y, L, M, J, diagIlincJ, dlamJ)

            if nargin < 12
                sysx = v;
                if nargin < 5
                    idx = []
                else
                    idx = z;
                end
                [Y, L, M, N, i] = obj.get_params(idx);
                [v, z, vi] = obj.x2vz(x, sysx, idx);

                %% compute linear current injections
                Ilin = Y*v + L*z + i;

                n  = length(v);     %% number of all port voltages
                ni = length(vi);    %% number of selected port voltages
                diagvi = sparse(1:ni, 1:ni, vi, ni, ni);
                if isempty(idx)     %% all ports
                    diagIlincJ = sparse(1:n, 1:n, conj(Ilin), n, n);
                    dlamJ      = sparse(1:n, 1:n, lam, n, n);
                    J          = sparse(1:n, 1:n, 1, n, n);
                else                %% selected ports
                    diagIlincJ = sparse(1:ni, idx, conj(Ilin), ni, n);
                    dlamJ      = sparse(1:ni, idx, lam, ni, n);
                    J          = sparse(1:ni, idx, 1, ni, n);
                end
            else
                n  = length(v);     %% number of all port voltages
                ni = length(lam);
            end

            %% intermediate terms
            diagv  = sparse(1:n, 1:n, v, n, n);

            A = diagvi * diagIlincJ;
            B = diagvi * conj(Y);
            C = B * conj(diagv);
            D = sparse(1:n, 1:n, 1 ./ abs(v), n, n);
            E = diagvi * conj(L);
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

        function [Svazr, Svazi, Svmzr, Svmzi] = port_inj_power_hess_vz(obj, x, lam, v, z, diagvi, L, dlamJ)
            % [Svazr, Svazi, Svmzr, Svmzi] = obj.port_inj_power_hess_vz(x, lam)
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
                [Y, L] = obj.get_params(idx);
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
            D = sparse(1:n, 1:n, 1 ./ abs(v), n, n);
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
