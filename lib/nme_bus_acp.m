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
            Vamin = -Inf(nb, 1);
            Vamax =  Inf(nb, 1);
            k = find(dme.isref);
            Vamin(k) = dme.Va0(k);
            Vamax(k) = dme.Va0(k);

            nm.add_var('va', 'Va', nb, dme.Va0, Vamin, Vamax);
            nm.add_var('vm', 'Vm', nb, dme.Vm0, dme.Vmin, dme.Vmax);
        end

        %%-----  PF methods  -----
        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            %% complex bus voltages
            nn = nm.get_idx('node');
            V = nm.soln.v(nn.i1.bus:nn.iN.bus);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.update(dm, 'Va', angle(V), 'Vm', abs(V));
        end

        %%-----  OPF methods  -----
        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
            %% complex bus voltages
            nn = nm.get_idx('node');
            V = nm.soln.v(nn.i1.bus:nn.iN.bus);

            %% shadow prices on voltage magnitudes
            [vv, nne] = mm.get_idx('var', 'nle');
            lambda = mm.soln.lambda;
            muVmin = lambda.lower(vv.i1.Vm:vv.iN.Vm);
            muVmax = lambda.upper(vv.i1.Vm:vv.iN.Vm);

            %% shadow prices on node power balance
            [lamP, lamQ] = nm.opf_node_power_balance_prices(mm);
            lamP = lamP(nn.i1.bus:nn.iN.bus);   %% for bus nodes only
            lamQ = lamQ(nn.i1.bus:nn.iN.bus);   %% for bus nodes only

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.update(dm, 'Va', angle(V), 'Vm', abs(V), ...
                'lamP', lamP, 'lamQ', lamQ, 'muVmin', muVmin, 'muVmax', muVmax);
        end
    end     %% methods
end         %% classdef
