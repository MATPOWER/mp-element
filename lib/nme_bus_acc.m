classdef nme_bus_acc < nme_bus & mp_form_acc

%   MATPOWER
%   Copyright (c) 2018-2022, Power Systems Engineering Research Center (PSERC)
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
    end     %% methods
end         %% classdef
