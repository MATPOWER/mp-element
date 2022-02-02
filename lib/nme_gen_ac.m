classdef nme_gen_ac < nme_gen% & mp_form_ac

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties
    
    methods
        function obj = add_zvars(obj, nm, dm, idx)
            ng = obj.nk;
            dme = obj.data_model_element(dm);
            nm.add_var('zr', 'Pg', ng, dme.pg_start, dme.pg_lb, dme.pg_ub);
            nm.add_var('zi', 'Qg', ng, dme.qg_start, dme.qg_lb, dme.qg_ub);
        end

        function obj = build_params(obj, nm, dm)
            build_params@nme_gen(obj, nm, dm);      %% call parent
            obj.N = -speye(obj.nk * obj.nz);
        end

        %%-----  PF methods  -----
        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            %% generator active power
            ss = nm.get_idx('state');
            Sg = nm.soln.z(ss.i1.gen:ss.iN.gen) * dm.base_mva;

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pg(dme.on) = real(Sg);
            dme.tab.qg(dme.on) = imag(Sg);
        end

        %%-----  OPF methods  -----
        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
            %% generator active power
            ss = nm.get_idx('state');
            Sg = nm.soln.z(ss.i1.gen:ss.iN.gen) * dm.base_mva;
            vm_setpoint = abs(obj.C' * nm.soln.v);

            %% shadow prices on generator limits
            [vv, ll] = mm.get_idx();
            lambda = mm.soln.lambda;
            mu_pg_ub = lambda.upper(vv.i1.Pg:vv.iN.Pg);
            mu_pg_lb = lambda.lower(vv.i1.Pg:vv.iN.Pg);
            mu_qg_ub = lambda.upper(vv.i1.Qg:vv.iN.Qg);
            mu_qg_lb = lambda.lower(vv.i1.Qg:vv.iN.Qg);

            %% gen PQ capability curve multipliers - based on update_mupq()
            if ll.N.PQh > 0 || ll.N.PQl > 0
                d = mm.get_userdata('Apqdata');

                %% combine original limit multipliers into single value
                muP = mu_pg_ub - mu_pg_lb;
                muQ = mu_qg_ub - mu_qg_lb;

                %% add P and Q components of multipliers on upper sloped constraint
                if ~isempty(d.ipqh)
                    mu_PQh = lambda.mu_l(ll.i1.PQh:ll.iN.PQh) - lambda.mu_u(ll.i1.PQh:ll.iN.PQh);
                    muP(d.ipqh) = muP(d.ipqh) - mu_PQh .* d.h(:,1);
                    muQ(d.ipqh) = muQ(d.ipqh) - mu_PQh .* d.h(:,2);
                end

                %% add P and Q components of multipliers on lower sloped constraint
                if ~isempty(d.ipql)
                    mu_PQl = lambda.mu_l(ll.i1.PQl:ll.iN.PQl) - lambda.mu_u(ll.i1.PQl:ll.iN.PQl);
                    muP(d.ipql) = muP(d.ipql) - mu_PQl .* d.l(:,1);
                    muQ(d.ipql) = muQ(d.ipql) - mu_PQl .* d.l(:,2);
                end

                %% split back into upper and lower multipliers based on sign
                mu_pg_ub = (muP > 0) .*  muP;
                mu_pg_lb = (muP < 0) .* -muP;
                mu_qg_ub = (muQ > 0) .*  muQ;
                mu_qg_lb = (muQ < 0) .* -muQ;
            end

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pg(dme.on) = real(Sg);
            dme.tab.qg(dme.on) = imag(Sg);
            dme.tab.vm_setpoint(dme.on) = vm_setpoint;
            dme.tab.mu_pg_lb(dme.on) = mu_pg_lb / dm.base_mva;
            dme.tab.mu_pg_ub(dme.on) = mu_pg_ub / dm.base_mva;
            dme.tab.mu_qg_lb(dme.on) = mu_qg_lb / dm.base_mva;
            dme.tab.mu_qg_ub(dme.on) = mu_qg_ub / dm.base_mva;
        end
    end     %% methods
end         %% classdef
