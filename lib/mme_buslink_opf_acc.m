classdef mme_buslink_opf_acc < mme_buslink

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'buslink';
%     end
    
    methods
        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            nme = obj.network_model_element(nm);

            %% voltage equality constraints
            [A, b_va, b_vm] = nme.voltage_constraints();
            vs = struct('name', {'Vr', 'Vr3', 'Vr3', 'Vr3', ...
                                 'Vi', 'Vi3', 'Vi3', 'Vi3'}, ...
                        'idx', {{}, {1}, {2}, {3}, {}, {1}, {2}, {3}});

            fcn_va = @(xx)opf_va_fcn(nme, xx, A, b_va);
            hess_va = @(xx, lam)opf_va_hess(nme, xx, lam, A);
            mm.add_nln_constraint('buslink_va', length(b_va), 1, fcn_va, hess_va, vs);

            fcn_vm = @(xx)opf_vm2_fcn(nme, xx, A, b_vm);
            hess_vm = @(xx, lam)opf_vm2_hess(nme, xx, lam, A);
            mm.add_nln_constraint('buslink_vm', length(b_vm), 1, fcn_vm, hess_vm, vs);
        end

        function x0 = opf_interior_x0(obj, mm, nm, dm, x0)
        end
    end     %% methods
end         %% classdef
