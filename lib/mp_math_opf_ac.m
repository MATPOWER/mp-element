classdef mp_math_opf_ac < mp_math_opf
%MP_MATH_OPF_AC  MATPOWER mathematical model for AC optimal power flow (OPF) problem.
%   ?
%
%   MP_MATH_OPF_AC ... power flow ...
%
%   Properties
%       ? - ?
%
%   Methods
%       ?

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
        function [g, dg] = nodal_current_balance_fcn(obj, nm, x_)
            if nargout > 1
                [G, Gv1, Gv2, Gzr, Gzi] = nm.nodal_complex_current_balance(x_);
                Gx = [Gv1 Gv2 Gzr Gzi];
                dg = [  real(Gx);       %% Re{I} mismatch w.r.t v1, v2, zr, zi
                        imag(Gx)    ];  %% Im{I} mismatch w.r.t v1, v2, zr, zi
            else
                G = nm.nodal_complex_current_balance(x_);
            end
            g = [ real(G);              %% real current mismatch
                  imag(G) ];            %% imaginary current mismatch
        end

        function [g, dg] = nodal_power_balance_fcn(obj, nm, x_)
            if nargout > 1
                [G, Gv1, Gv2, Gzr, Gzi] = nm.nodal_complex_power_balance(x_);
                Gx = [Gv1 Gv2 Gzr Gzi];
                dg = [  real(Gx);       %% P mismatch w.r.t v1, v2, zr, zi
                        imag(Gx)    ];  %% Q mismatch w.r.t v1, v2, zr, zi
            else
                G = nm.nodal_complex_power_balance(x_);
            end
            g = [ real(G);              %% active power (P) mismatch
                  imag(G) ];            %% reactive power (Q) mismatch
        end

        function d2G = nodal_current_balance_hess(obj, nm, x_, lam)
            nlam = length(lam) / 2;
            lamIr = lam(1:nlam);
            lamIi = lam((1:nlam)+nlam);

            d2Gr = nm.nodal_complex_current_balance_hess(x_, lamIr);
            d2Gi = nm.nodal_complex_current_balance_hess(x_, lamIi);

            d2G = real(d2Gr) + imag(d2Gi);
        end

        function d2G = nodal_power_balance_hess(obj, nm, x_, lam)
            nlam = length(lam) / 2;
            lam_p = lam(1:nlam);
            lam_q = lam((1:nlam)+nlam);

            d2Gr = nm.nodal_complex_power_balance_hess(x_, lam_p);
            d2Gi = nm.nodal_complex_power_balance_hess(x_, lam_q);

            d2G = real(d2Gr) + imag(d2Gi);
        end

        function add_system_costs(obj, nm, dm, mpopt)
            %% legacy user-defined costs
            if isfield(dm.userdata, 'legacy_opf_user_mods')
                obj.add_legacy_user_costs(nm, dm, 0);
            end
        end

        function add_legacy_user_constraints(obj, nm, dm, mpopt)
            %% call parent
            add_legacy_user_constraints@mp_math_opf(obj, nm, dm, mpopt);

            if ~isempty(dm.userdata.legacy_opf_user_mods)
                uc = dm.userdata.legacy_opf_user_mods.nlc;
                for k = 1:length(uc)
                    obj.add_nln_constraint(uc{k}{:});
                end
            end
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            opt = mpopt2nlpopt(mpopt, obj.problem_type());

            if mpopt.opf.start < 2
                %% initialize interior point
                opt.x0 = nm.opf_interior_x0(obj, nm, dm);
            end
        end
    end     %% methods
end         %% classdef