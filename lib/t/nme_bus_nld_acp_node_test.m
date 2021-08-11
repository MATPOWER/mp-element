classdef nme_bus_nld_acp_node_test < nme_bus_acp

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
        %% constructor
        function obj = nme_bus_nld_acp_node_test()
            obj@nme_bus_acp();
            obj.name = 'bus_nld';
        end

        function obj = add_vvars(obj, nm, dm, idx)
            dme = obj.data_model_element(dm);
            nb = obj.nk;

            %% prepare angle bounds for ref buses
            Vamin = -Inf(nb, 1);
            Vamax =  Inf(nb, 1);
            k = find(dme.type == NODE_TYPE.REF);
            Vamin(k) = dme.Va0(k);
            Vamax(k) = dme.Va0(k);

            nm.add_var('va', ['va_' obj.name], nb, dme.Va0, Vamin, Vamax);
            nm.add_var('vm', ['vm_' obj.name], nb, dme.Vm0, dme.Vmin, dme.Vmax);
        end

        %%-----  PF methods  -----
        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            %% complex bus voltages
            nn = nm.get_idx('node');
            V = nm.soln.v(nn.i1.(obj.name):nn.iN.(obj.name));

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.update(dm, 'Va', angle(V), 'Vm', abs(V));
        end

        %%-----  OPF methods  -----
        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
            %% complex bus voltages
            nn = nm.get_idx('node');
            V = nm.soln.v(nn.i1.(obj.name):nn.iN.(obj.name));

            %% shadow prices on voltage magnitudes
            vv = mm.get_idx('var');
            lambda = mm.soln.lambda;
            vVm = ['vm_' obj.name];
            muVmin = lambda.lower(vv.i1.(vVm):vv.iN.(vVm));
            muVmax = lambda.upper(vv.i1.(vVm):vv.iN.(vVm));

            %% shadow prices on node power balance
            [lamP, lamQ] = nm.opf_node_power_balance_prices(mm);
            lamP = lamP(nn.i1.(obj.name):nn.iN.(obj.name)); %% for (obj.name) nodes only
            lamQ = lamQ(nn.i1.(obj.name):nn.iN.(obj.name)); %% for (obj.name) nodes only

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.update(dm, 'Va', angle(V), 'Vm', abs(V), ...
                'lamP', lamP, 'lamQ', lamQ, 'muVmin', muVmin, 'muVmax', muVmax);
        end
    end     %% methods
end         %% classdef
