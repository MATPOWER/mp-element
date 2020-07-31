classdef acc_bus < mp_bus & mp_model_acc

%   MATPOWER
%   Copyright (c) 2018-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Baljinnyam Sereeter, Delft University of Technology
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'bus';
%     end
    
    methods
        function obj = add_vvars(obj, nm, mpc, idx)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            nb = obj.nk;
            v0 = mpc.bus(:, VM) .* exp(1j * mpc.bus(:, VA) * pi/180);
            vclim = 1.1 * mpc.bus(:, VMAX);
            nm.add_var('vr', 'Vr', nb, real(v0), -vclim, vclim);
            nm.add_var('vi', 'Vi', nb, imag(v0), -vclim, vclim);
        end

        function [g, dg] = va_fcn(obj, xx, idx, lim)
            %% lim can be a scalar value for equality constraint or
            %% upper bound on va, or a cell array with {vamin, vamax}
            %% for double bounds

            %% unpack data
            [vr, vi] = deal(xx{:});

            %% compute voltage angle mismatch
            va = angle(vr(idx) + 1j* vi(idx));
            if iscell(lim)
                g = [ lim{1} - va;
                      va - lim{2}   ];
            else
                g = va - lim;
            end

            if nargout > 1
                %% compute partials of voltage angle w.r.t vr and vi
                n = length(idx);
                nn = length(vr);
                vm2 = vr(idx).^2 + vi(idx).^2;
                dva_dvr = sparse(1:n, idx, -vi(idx) ./ vm2, n, nn);
                dva_dvi = sparse(1:n, idx,  vr(idx) ./ vm2, n, nn);
                if iscell(lim)
                    dg = [ -dva_dvr -dva_dvi;       %% va w.r.t vr, vi
                            dva_dvr  dva_dvi  ];    %% va w.r.t vr, vi
                else
                    dg = [dva_dvr dva_dvi];         %% va w.r.t vr, vi
                end
            end
        end

        function d2G = va_hess(obj, xx, lam, idx)
            %% unpack data
            [vr, vi] = deal(xx{:});
            n = length(idx);
            nn = length(vr);

            %% evaluate Hessian of voltage angle function
            vvr = vr(idx);
            vvi = vi(idx);
            vvr2 = vvr.^2;
            vvi2 = vvi.^2;
            vvm4 = (vvr2 + vvi2).^2;
            if length(lam) == n     %% upper bound or equality
                lamvm4 = lam ./ vvm4;
            else                    %% doubly bounded (use lam_ub-lam_lb)
                lamvm4 = (lam(n+1:2*n) - lam(1:n)) ./ vvm4;
            end
            d2vref_rr = sparse(idx, idx, 2 * lamvm4 .*  vvr .* vvi,   nn, nn);
            d2vref_ri = sparse(idx, idx,     lamvm4 .* (vvi2 - vvr2), nn, nn);

            %% construct Hessian
            d2G = [ d2vref_rr  d2vref_ri;
                    d2vref_ri -d2vref_rr ];
        end

        function [g, dg] = vm2_fcn(obj, xx, idx, lim)
            %% lim can be a scalar value for equality constraint or
            %% upper bound on vm^2, or a cell array with {vm2min, vm2max}
            %% for double bounds

            %% unpack data
            [vr, vi] = deal(xx{:});

            %% compute voltage magnitude^2 mismatch
            vm2 = vr(idx).^2 + vi(idx).^2;
            if iscell(lim)
                g = [ lim{1} - vm2;
                      vm2 - lim{2}   ];
            else
                g = vm2 - lim;
            end

            if nargout > 1
                %% compute partials of voltage magnitude^2 w.r.t vr and vi
                n = length(idx);
                nn = length(vr);
                dvm_dvr = sparse(1:n, idx, 2 * vr(idx), n, nn);
                dvm_dvi = sparse(1:n, idx, 2 * vi(idx), n, nn);
                if iscell(lim)
                    dg = [ -dvm_dvr -dvm_dvi;
                            dvm_dvr  dvm_dvi  ];    %% vm2 w.r.t vr, vi
                else
                    dg = [ dvm_dvr dvm_dvi ];       %% vm2 w.r.t vr, vi
                end
            end
        end

        function d2G = vm2_hess(obj, xx, lam, idx)
            %% unpack data
            [vr, vi] = deal(xx{:});
            n = length(idx);
            nn = length(vr);
            
            %% evaluate Hessian of voltage magnitude^2 function
            if length(lam) == n     %% upper bound or equality
                dlam = sparse(idx, idx, 2*lam, nn, nn);
            else                    %% doubly bounded (use lam_ub-lam_lb)
                dlam = sparse(idx, idx, 2*(lam(n+1:2*n) - lam(1:n)), nn, nn);
            end

            %% construct Hessian
            zz = sparse(nn, nn);
            d2G = [dlam zz; zz dlam];
        end

        function add_opf_constraints(obj, nm, om, mpc, mpopt)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            %% voltage angle reference constraint
            refs = find(mpc.bus(:, BUS_TYPE) == REF);
            varef = mpc.bus(refs, VA) * pi/180;
            fcn_vref = @(xx)va_fcn(obj, xx, refs, varef);
            hess_vref = @(xx, lam)va_hess(obj, xx, lam, refs);
            om.add_nln_constraint('Vref', length(refs), 1, fcn_vref, hess_vref, {'Vr', 'Vi'});

            %% fixed voltage magnitudes
            veq = find(mpc.bus(:, VMIN) == mpc.bus(:, VMAX));
            nveq = length(veq);
            if nveq
                fcn_vm2eq = @(xx)vm2_fcn(obj, xx, veq, mpc.bus(veq, VMAX).^2);
                hess_vm2eq = @(xx, lam)vm2_hess(obj, xx, lam, veq);
                om.add_nln_constraint('Veq', nveq, 1, fcn_vm2eq, hess_vm2eq, {'Vr', 'Vi'});
            end
            om.userdata.veq = veq;

            %% voltage magnitude limits
            viq = find(mpc.bus(:, VMIN) ~= mpc.bus(:, VMAX));
            nviq = length(viq);
            if nviq
                fcn_vlim = @(xx)vm2_fcn(obj, xx, viq, ...
                        {mpc.bus(viq, VMIN).^2, mpc.bus(viq, VMAX).^2} );
                hess_vlim = @(xx, lam)vm2_hess(obj, xx, lam, viq);
                om.add_nln_constraint({'Vmin', 'Vmax'}, [nviq;nviq], 0, fcn_vlim, hess_vlim, {'Vr', 'Vi'});
            end
            om.userdata.viq = viq;
        end
    end     %% methods
end         %% classdef
