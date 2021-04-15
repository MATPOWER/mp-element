classdef nme_bus_dc < nme_bus & mp_form_dc

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
        end

        %%-----  PF methods  -----
        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            %% bus voltage angles
            nn = nm.get_idx('node');
            Va = nm.soln.v(nn.i1.bus:nn.iN.bus);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.update(dm, 'Va', Va, 'Vm', 1);
        end

        %%-----  OPF methods  -----
        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
            %% bus voltage angles
            nn = nm.get_idx('node');
            Va = nm.soln.v(nn.i1.bus:nn.iN.bus);

            %% shadow prices on node power balance
            ll = mm.get_idx('lin');
            lambda = mm.soln.lambda;
            lamP =  lambda.mu_u(ll.i1.Pmis:ll.iN.Pmis) - ...
                    lambda.mu_l(ll.i1.Pmis:ll.iN.Pmis);
            lamP = lamP(nn.i1.bus:nn.iN.bus);   %% for bus nodes only

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.update(dm, 'Va', Va, 'Vm', 1, 'lamP', lamP);
        end
    end     %% methods
end         %% classdef
