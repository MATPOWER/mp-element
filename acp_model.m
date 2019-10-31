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

        function [Iva, Ivm, Izr, Izi] = port_inj_current_jac(obj, ...
                n, vv, Y, M, diagv, invdiagvic, diagSlinc)
            % [Iva, Ivm] = obj.port_inj_current_jac(...)
            % [Iva, Ivm, Izr, Izi] = obj.port_inj_current_jac(...)

            %% intermediate terms
            A = invdiagvic * (diagSlinc - conj(M * diagv));
            D = sparse(1:n, 1:n, 1 ./ abs(vv), n, n);

            %% linear current term
            Ilin_vm = Y * diagv;
            Ilin_va = 1j * Ilin_vm;
            Ilin_vm = Ilin_vm * D;

            %% current from linear power term
            IS_va = 1j * A;
            IS_vm = -A * D;

            %% combine
            Iva = Ilin_va + IS_va;
            Ivm = Ilin_vm + IS_vm;
        end

        function [Sva, Svm] = port_inj_power_jac(obj, ...
                n, vv, Y, M, diagv, diagvi, diagIlinc)
            % [Sva, Svm] = obj.port_inj_power_jac(...)

            %% intermediate terms
            A = diagvi * diagIlinc;
            B = diagvi * conj(Y);
            C = B * conj(diagv);
            D = sparse(1:n, 1:n, 1 ./ abs(vv), n, n);

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
    end     %% methods
end         %% classdef
