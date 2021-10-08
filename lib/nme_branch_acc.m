classdef nme_branch_acc < nme_branch_ac & mp_form_acc

%   MATPOWER
%   Copyright (c) 2018-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.


    methods
        function [h, dh] = ang_diff_fcn(obj, xx, Aang, lang, uang)
            if nargout > 1
                [va, dva] = obj.va_fcn(xx, [], 0);
                Adva = Aang * dva;
                dh = [ -Adva; Adva ];
            else
                va = obj.va_fcn(xx, [], 0);
            end
            Ax = Aang * va;
            h = [ lang - Ax;
                  Ax - uang ];
        end

        function d2H = ang_diff_hess(obj, xx, lambda, Aang)
            %% unpack data
            [vr, vi] = deal(xx{:});
            nn = length(vr);

            %% evaluate Hessian of voltage magnitude^2 function
            n = length(lambda) / 2;
            if n
                lam = lambda(n+1:2*n) - lambda(1:n);    %% lam_ub - lam_lb
            else
                lam = zeros(0,1);
            end

            d2H = obj.va_hess(xx, Aang' * lam, []);
        end

        function [mu_vad_lb, mu_vad_ub] = opf_branch_ang_diff_prices(obj, mm)
            %% shadow prices on angle difference limits
            iang = mm.userdata.ang_diff_constrained_branch_idx;
            mu_vad_lb = zeros(obj.nk, 1);
            mu_vad_ub = mu_vad_lb;
            if length(iang)
                nni = mm.get_idx('nli');
                lambda = mm.soln.lambda;
                mu_vad_ub(iang) = lambda.ineqnonlin(nni.i1.angU:nni.iN.angU);
                mu_vad_lb(iang) = lambda.ineqnonlin(nni.i1.angL:nni.iN.angL);
            end
        end
    end     %% methods
end         %% classdef
