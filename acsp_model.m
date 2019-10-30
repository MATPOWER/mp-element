classdef acsp_model < ac_model
%ACSP_MODEL  MATPOWER Model class for AC-power-polar models.
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
        function name = model_name(obj)
            name = 'AC-power-polar model';
        end
        function tag = model_tag(obj)
            tag = 'acsp';
        end
        function vtypes = model_vvars(obj)
            vtypes = {'va', 'vm'};
        end

        function [Sva, Svm, Szr, Szi] = port_inj_power_jac(obj, ...
                n, ni, vv, Y, L, M, N, diagv, diagvi, diagIlinc)
            % [Sva, Svm] = obj.port_inj_power(...)
            % [Sva, Svm, Szr, Szi] = obj.port_inj_power(...)

            %% intermediate terms
            A = diagvi * diagIlinc;
            B = diagvi * conj(Y);
            C = B * conj(diagv);
            D = sparse(1:n, 1:n, 1 ./ abs(vv), n, n);

            %% linear power term
            Slin_va = 1j * M * diagv;
            Slin_vm =  M * D * diagv;

            %% power from linear current term
            SI_va = 1j * (A - C);
            SI_vm = (A + C) * D;

            %% combine
            Sva = Slin_va + SI_va;
            Svm = Slin_vm + SI_vm;

            if nargout > 2
                %% intermediate terms
                E = diagvi * conj(L);

                %% linear power term
                Slin_zr = N;
                Slin_zi = 1j * N;

                %% power from linear current term
                SI_zr = E;
                SI_zi = 1j * E;

                %% combine
                Szr = Slin_zr + SI_zr;
                Szi = Slin_zi + SI_zi;
            end
        end

    end     %% methods
end         %% classdef
