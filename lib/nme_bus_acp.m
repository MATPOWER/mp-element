classdef nme_bus_acp < nme_bus & mp_form_acp

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
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

            %% prepare angle bounds for ref buses
            va_lb = -Inf(nb, 1);
            va_ub =  Inf(nb, 1);
            k = find(dme.type == NODE_TYPE.REF);
            va_lb(k) = dme.va_start(k);
            va_ub(k) = dme.va_start(k);

            nm.add_var('va', 'Va', nb, dme.va_start, va_lb, va_ub);
            nm.add_var('vm', 'Vm', nb, dme.vm_start, dme.vm_lb, dme.vm_ub);
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
            vv = mm.get_idx('var');
            lambda = mm.soln.lambda;
            mu_vm_lb = lambda.lower(vv.i1.Vm:vv.iN.Vm);
            mu_vm_ub = lambda.upper(vv.i1.Vm:vv.iN.Vm);

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
            varef1 = nm.opf_interior_va(mm, dm);
            vm = obj.opf_interior_vm(mm, nm, dm);
            x0(vv.i1.Va:vv.iN.Va) = varef1; %% angles set to 1st ref angle
            x0(vv.i1.Vm:vv.iN.Vm) = vm;     %% voltage magnitudes
        end
    end     %% methods
end         %% classdef
