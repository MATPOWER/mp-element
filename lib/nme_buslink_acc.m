classdef nme_buslink_acc < nme_buslink & mp_form_acc

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        %%-----  PF methods  -----
        function [g, dg] = pf_va_fcn(obj, xx, A, b)
            %% unpack data
            vr = vertcat(xx{1:8});
            vi = vertcat(xx{9:16});

            if nargout > 1
                [va, dva] = obj.va_fcn({vr, vi}, [], 0);
                dg = A * dva;
            else
                va = obj.va_fcn({vr, vi}, [], 0);
            end
            g = A * va - b;
        end

        function [g, dg] = pf_vm_fcn(obj, xx, A, b)
            %% unpack data
            vr = vertcat(xx{1:4});
            vi = vertcat(xx{5:8});

            if nargout > 1
                n = length(vr);
                [vm, dvm] = obj.vm2_fcn({vr, vi}, [], 0);
                vm_inv = 1./vm;
                dg = 0.5 * A * dvm * spdiags([vm_inv; vm_inv], 0, 2*n, 2*n);
            else
                vm = obj.vm2_fcn({vr, vi}, [], 0);
            end
            g = A * sqrt(vm) - b;
        end

        function obj = pf_add_constraints(obj, mm, nm, dm, mpopt)
            %% add constraints for matching
            %%  voltage angles at pv and pq nodes
            %%  voltage magnitudes at pq nodes
            [A_va_pq, A_va_pv, b_va, A_vm, b_vm] = pf_voltage_constraints(obj, mm.aux_data);

            %% prep variable set structs
            vs_va = struct('name', {'Vr_pq', 'Vr3_pq', 'Vr3_pq', 'Vr3_pq', ...
                                    'Vr_pv', 'Vr3_pv', 'Vr3_pv', 'Vr3_pv', ...
                                    'Vi_pq', 'Vi3_pq', 'Vi3_pq', 'Vi3_pq', ...
                                    'Vi_pv', 'Vi3_pv', 'Vi3_pv', 'Vi3_pv'}, ...
                            'idx', {{}, {1}, {2}, {3}, {}, {1}, {2}, {3}, ...
                                    {}, {1}, {2}, {3}, {}, {1}, {2}, {3}} );
            vs_vm = struct('name', {'Vr_pq', 'Vr3_pq', 'Vr3_pq', 'Vr3_pq', ...
                                    'Vi_pq', 'Vi3_pq', 'Vi3_pq', 'Vi3_pq'}, ...
                            'idx', {{}, {1}, {2}, {3}, {}, {1}, {2}, {3}} );

            fcn_va = @(xx)pf_va_fcn(obj, xx, [A_va_pq A_va_pv], b_va);
            mm.add_nln_constraint('buslink_va', length(b_va), 1, fcn_va, [], vs_va);

            fcn_vm = @(xx)pf_vm_fcn(obj, xx, A_vm, b_vm);
            mm.add_nln_constraint('buslink_vm', length(b_vm), 1, fcn_vm, [], vs_vm);

            %% call parent
            pf_add_constraints@nme_buslink(obj, mm, nm, dm, mpopt);
        end

        %%-----  OPF methods  -----
        function [g, dg] = opf_va_fcn(obj, xx, A, b)
            %% unpack data
            vr = vertcat(xx{1:4});
            vi = vertcat(xx{5:8});

            if nargout > 1
                [va, dva] = obj.va_fcn({vr, vi}, [], 0);
                dg = A * dva;
            else
                va = obj.va_fcn({vr, vi}, [], 0);
            end
            g = A * va - b;
        end

        function d2G = opf_va_hess(obj, xx, lam, A)
            %% unpack data
            vr = vertcat(xx{1:4});
            vi = vertcat(xx{5:8});

            d2G = obj.va_hess({vr, vi}, A' * lam, []);
        end

        function [g, dg] = opf_vm2_fcn(obj, xx, A, b)
            %% unpack data
            vr = vertcat(xx{1:4});
            vi = vertcat(xx{5:8});

            if nargout > 1
                [vm, dvm] = obj.vm2_fcn({vr, vi}, [], 0);
                dg = A * dvm;
            else
                vm = obj.vm2_fcn({vr, vi}, [], 0);
            end
            g = A * vm - b;
        end

        function d2G = opf_vm2_hess(obj, xx, lam, A)
            %% unpack data
            vr = vertcat(xx{1:4});
            vi = vertcat(xx{5:8});

            d2G = obj.vm2_hess({vr, vi}, A' * lam, []);
        end

        function obj = opf_add_constraints(obj, mm, nm, dm, mpopt)
            %% voltage equality constraints
            [A, b_va, b_vm] = obj.voltage_constraints();
            vs = struct('name', {'Vr', 'Vr3', 'Vr3', 'Vr3', ...
                                 'Vi', 'Vi3', 'Vi3', 'Vi3'}, ...
                        'idx', {{}, {1}, {2}, {3}, {}, {1}, {2}, {3}});

            fcn_va = @(xx)opf_va_fcn(obj, xx, A, b_va);
            hess_va = @(xx, lam)opf_va_hess(obj, xx, lam, A);
            mm.add_nln_constraint('buslink_va', length(b_va), 1, fcn_va, hess_va, vs);

            fcn_vm = @(xx)opf_vm2_fcn(obj, xx, A, b_vm);
            hess_vm = @(xx, lam)opf_vm2_hess(obj, xx, lam, A);
            mm.add_nln_constraint('buslink_vm', length(b_vm), 1, fcn_vm, hess_vm, vs);
        end
    end     %% methods
end         %% classdef
