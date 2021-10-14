classdef mme_bus_opf_acc < mme_bus

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'bus';
%     end
    
    methods
        function add_constraints(obj, mm, nm, dm, mpopt)
            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            %% voltage angle reference constraint
            ref = find(nme.node_types(nm, dm) == NODE_TYPE.REF);
            varef = dme.va_start(ref);
            fcn_vref = @(xx)va_fcn(nme, xx, ref, varef);
            hess_vref = @(xx, lam)va_hess(nme, xx, lam, ref);
            mm.add_nln_constraint('Vref', length(ref), 1, fcn_vref, hess_vref, {'Vr', 'Vi'});

            %% fixed voltage magnitudes
            veq = find(dme.vm_lb == dme.vm_ub);
            nveq = length(veq);
            if nveq
                fcn_vm2eq = @(xx)vm2_fcn(nme, xx, veq, dme.vm_ub(veq).^2);
                hess_vm2eq = @(xx, lam)vm2_hess(nme, xx, lam, veq);
                mm.add_nln_constraint('Veq', nveq, 1, fcn_vm2eq, hess_vm2eq, {'Vr', 'Vi'});
            end
            mm.userdata.veq = veq;

            %% voltage magnitude limits
            viq = find(dme.vm_lb ~= dme.vm_ub);
            nviq = length(viq);
            if nviq
                fcn_vlim = @(xx)vm2_fcn(nme, xx, viq, ...
                        {dme.vm_lb(viq).^2, dme.vm_ub(viq).^2} );
                hess_vlim = @(xx, lam)vm2_hess(nme, xx, lam, viq);
                mm.add_nln_constraint({'Vmin', 'Vmax'}, [nviq;nviq], 0, fcn_vlim, hess_vlim, {'Vr', 'Vi'});
            end
            mm.userdata.viq = viq;
        end
    end     %% methods
end         %% classdef
