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
            V0 = dme.vm_start .* exp(1j * dme.va_start);
            vclim = 1.1 * dme.vm_ub;

            nm.add_var('vr', 'Vr', nb, real(V0), -vclim, vclim);
            nm.add_var('vi', 'Vi', nb, imag(V0), -vclim, vclim);
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
        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
            %% complex bus voltages
            nn = nm.get_idx('node');
            V = nm.soln.v(nn.i1.bus:nn.iN.bus);

            %% shadow prices on voltage magnitudes
            [nne, nni] = mm.get_idx('nle', 'nli');
            lambda = mm.soln.lambda;
            mu_vm_lb = zeros(nn.N.bus, 1);  %% init to all 0
            mu_vm_ub = mu_vm_lb;            %% init to all 0
            if mm.userdata.veq
                lam = lambda.eqnonlin(nne.i1.Veq:nne.iN.Veq);
                lam_p = zeros(size(lam));
                lam_n = zeros(size(lam));
                lam_p(lam > 0) =  lam(lam > 0);
                lam_n(lam < 0) = -lam(lam < 0);
                mu_vm_lb(mm.userdata.veq) = lam_n;
                mu_vm_ub(mm.userdata.veq) = lam_p;
            end
            mu_vm_lb(mm.userdata.viq) = lambda.ineqnonlin(nni.i1.Vmin:nni.iN.Vmin);
            mu_vm_ub(mm.userdata.viq) = lambda.ineqnonlin(nni.i1.Vmax:nni.iN.Vmax);

            vm = abs(V);
            mu_vm_lb = mu_vm_lb .* vm * 2;
            mu_vm_ub = mu_vm_ub .* vm * 2;

            %% shadow prices on node power balance
            [lam_p, lam_q] = nm.opf_node_power_balance_prices(mm);
            lam_p = lam_p(nn.i1.bus:nn.iN.bus);     %% for bus nodes only
            lam_q = lam_q(nn.i1.bus:nn.iN.bus);     %% for bus nodes only

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.va(dme.on) = angle(V) * 180/pi;
            dme.tab.vm(dme.on) = abs(V);
            dme.tab.lam_p(dme.on) = lam_p / dm.base_mva;
            dme.tab.lam_q(dme.on) = lam_q / dm.base_mva;
            dme.tab.mu_vm_lb(dme.on) = mu_vm_lb;
            dme.tab.mu_vm_ub(dme.on) = mu_vm_ub;
        end

        function x0 = opf_interior_x0(obj, mm, nm, dm, x0)
            vv = mm.get_idx();
            varef1 = mm.opf_interior_va(nm, dm);
            vm = obj.opf_interior_vm(mm, nm, dm);
            v_ = vm * exp(1j*varef1);
            x0(vv.i1.Vr:vv.iN.Vr) = real(v_);
            x0(vv.i1.Vi:vv.iN.Vi) = imag(v_);
        end
    end     %% methods
end         %% classdef
