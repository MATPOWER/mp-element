classdef nme_gen_dc < nme_gen & mp_form_dc

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'gen';
%     end
    
    methods
        function obj = add_zvars(obj, nm, dm, idx)
            ng = obj.nk;
            dme = obj.data_model_element(dm);
            nm.add_var('z', 'Pg', ng, dme.pg_start, dme.pg_lb, dme.pg_ub);
        end

        function obj = build_params(obj, nm, dm)
            build_params@nme_gen(obj, nm, dm);      %% call parent
            obj.K = -speye(obj.nk * obj.nz);
        end

        %%-----  PF methods  -----
        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            %% generator active power
            ss = nm.get_idx('state');
            pg = nm.soln.z(ss.i1.gen:ss.iN.gen) * dm.base_mva;

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pg(dme.on) = pg;
        end

        %%-----  OPF methods  -----
        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
            %% generator active power
            ss = nm.get_idx('state');
            pg = nm.soln.z(ss.i1.gen:ss.iN.gen) * dm.base_mva;

            %% shadow prices on generator limits
            vv = mm.get_idx();
            lambda = mm.soln.lambda;
            mu_pg_ub = lambda.upper(vv.i1.Pg:vv.iN.Pg);
            mu_pg_lb = lambda.lower(vv.i1.Pg:vv.iN.Pg);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pg(dme.on) = pg;
            dme.tab.mu_pg_lb(dme.on) = mu_pg_lb / dm.base_mva;
            dme.tab.mu_pg_ub(dme.on) = mu_pg_ub / dm.base_mva;
        end
    end     %% methods
end         %% classdef
