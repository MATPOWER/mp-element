classdef nme_bus_acc < nme_bus & mp_form_acc

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
        function obj = add_vvars(obj, nm, dm, idx)
            dme = obj.data_model_element(dm);
            nb = obj.nk;

            %% prepare box bounds for voltage coordinates
            v0 = dme.Vm0 .* exp(1j * dme.Va0);
            vclim = 1.1 * dme.Vmax;

            nm.add_var('vr', 'Vr', nb, real(v0), -vclim, vclim);
            nm.add_var('vi', 'Vi', nb, imag(v0), -vclim, vclim);
        end

        %%-----  PF methods  -----
        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            %% complex bus voltages
            nn = nm.get_idx('node');
            V = nm.soln.v(nn.i1.bus:nn.iN.bus);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.va(dme.on) = angle(V) * 180/pi;
            dme.tab.vm(dme.on) = abs(V);
        end

        %%-----  OPF methods  -----
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

        function opf_add_constraints(obj, mm, nm, dm, mpopt)
            %% voltage angle reference constraint
            dme = obj.data_model_element(dm);
            ref = find(obj.node_types(nm, dm) == NODE_TYPE.REF);
            varef = dme.Va0(ref);
            fcn_vref = @(xx)va_fcn(obj, xx, ref, varef);
            hess_vref = @(xx, lam)va_hess(obj, xx, lam, ref);
            mm.add_nln_constraint('Vref', length(ref), 1, fcn_vref, hess_vref, {'Vr', 'Vi'});

            %% fixed voltage magnitudes
            veq = find(dme.Vmin == dme.Vmax);
            nveq = length(veq);
            if nveq
                fcn_vm2eq = @(xx)vm2_fcn(obj, xx, veq, dme.Vmax(veq).^2);
                hess_vm2eq = @(xx, lam)vm2_hess(obj, xx, lam, veq);
                mm.add_nln_constraint('Veq', nveq, 1, fcn_vm2eq, hess_vm2eq, {'Vr', 'Vi'});
            end
            mm.userdata.veq = veq;

            %% voltage magnitude limits
            viq = find(dme.Vmin ~= dme.Vmax);
            nviq = length(viq);
            if nviq
                fcn_vlim = @(xx)vm2_fcn(obj, xx, viq, ...
                        {dme.Vmin(viq).^2, dme.Vmax(viq).^2} );
                hess_vlim = @(xx, lam)vm2_hess(obj, xx, lam, viq);
                mm.add_nln_constraint({'Vmin', 'Vmax'}, [nviq;nviq], 0, fcn_vlim, hess_vlim, {'Vr', 'Vi'});
            end
            mm.userdata.viq = viq;
        end

        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
            %% complex bus voltages
            nn = nm.get_idx('node');
            V = nm.soln.v(nn.i1.bus:nn.iN.bus);

            %% shadow prices on voltage magnitudes
            [nne, nni] = mm.get_idx('nle', 'nli');
            lambda = mm.soln.lambda;
            muVmin = zeros(nn.N.bus, 1);    %% init to all 0
            muVmax = muVmin;                %% init to all 0
            if mm.userdata.veq
                lam = lambda.eqnonlin(nne.i1.Veq:nne.iN.Veq);
                lam_p = zeros(size(lam));
                lam_n = zeros(size(lam));
                lam_p(lam > 0) =  lam(lam > 0);
                lam_n(lam < 0) = -lam(lam < 0);
                muVmin(mm.userdata.veq) = lam_n;
                muVmax(mm.userdata.veq) = lam_p;
            end
            muVmin(mm.userdata.viq) = lambda.ineqnonlin(nni.i1.Vmin:nni.iN.Vmin);
            muVmax(mm.userdata.viq) = lambda.ineqnonlin(nni.i1.Vmax:nni.iN.Vmax);

            Vm = abs(V);
            muVmin = muVmin .* Vm * 2;
            muVmax = muVmax .* Vm * 2;

            %% shadow prices on node power balance
            [lamP, lamQ] = nm.opf_node_power_balance_prices(mm);
            lamP = lamP(nn.i1.bus:nn.iN.bus);   %% for bus nodes only
            lamQ = lamQ(nn.i1.bus:nn.iN.bus);   %% for bus nodes only

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.va(dme.on) = angle(V) * 180/pi;
            dme.tab.vm(dme.on) = abs(V);
            dme.tab.lam_p(dme.on) = lamP / dm.baseMVA;
            dme.tab.lam_q(dme.on) = lamQ / dm.baseMVA;
            dme.tab.mu_vm_lb(dme.on) = muVmin;
            dme.tab.mu_vm_ub(dme.on) = muVmax;
        end
    end     %% methods
end         %% classdef
