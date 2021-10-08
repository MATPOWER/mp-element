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
    end     %% methods
end         %% classdef
