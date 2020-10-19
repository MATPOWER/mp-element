classdef mpe_branch_acc < mpe_branch_ac & mp_model_acc

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
            %% unpack data
            [vr, vi] = deal(xx{:});

            %% compute branch angle difference
            va = angle(vr + 1j* vi);
            Ax = Aang * va;
            h = [ lang - Ax;
                  Ax - uang ];

            if nargout > 1
                %% compute partials of branch angle difference w.r.t vr and vi
                nn = length(vr);
                vm2 = vr.^2 + vi.^2;
                Aangdva_dvr = Aang * sparse(1:nn, 1:nn, -vi./vm2, nn, nn);
                Aangdva_dvi = Aang * sparse(1:nn, 1:nn,  vr./vm2, nn, nn);
                dh = [ -Aangdva_dvr -Aangdva_dvi;   %% h w.r.t vr, vi
                        Aangdva_dvr  Aangdva_dvi ];
            end
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

            vr2 = vr.^2;
            vi2 = vi.^2;

            lam_vm4 = (Aang' * lam) ./ (vr2 + vi2).^2;

            h_rr = sparse(1:nn, 1:nn, 2 * lam_vm4 .*  vr .* vi,   nn, nn);
            h_ri = sparse(1:nn, 1:nn,     lam_vm4 .* (vi2 - vr2), nn, nn);

            %% construct Hessian
            d2H = [ h_rr  h_ri;
                    h_ri -h_rr ];
        end

        function add_opf_constraints(obj, nm, om, mpc, mpopt)
            %% call parent
            add_opf_constraints@mpe_branch_ac(obj, nm, om, mpc, mpopt);

            %% branch angle difference limits
            nb = size(mpc.bus, 1);
            [Aang, lang, uang, iang]  = makeAang(mpc.baseMVA, mpc.branch, nb, mpopt);
            nang = length(iang);
            if nang
                fcn_ang = @(xx)ang_diff_fcn(obj, xx, Aang, lang, uang);
                hess_ang = @(xx, lam)ang_diff_hess(obj, xx, lam, Aang);
                om.add_nln_constraint({'angL', 'angU'}, [nang;nang], 0, fcn_ang, hess_ang, {'Vr', 'Vi'});
            end
            om.userdata.iang = iang;
        end
    end     %% methods
end         %% classdef
