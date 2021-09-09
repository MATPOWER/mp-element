classdef nme_branch_dc < nme_branch & mp_form_dc

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'branch';
%     end
    
    methods
        function obj = build_params(obj, nm, dm)
            build_params@nme_branch(obj, nm, dm);   %% call parent

            dme = obj.data_model_element(dm);
            nl = obj.nk;

            tm = ones(nl, 1);           %% default tap ratio = 1
            i = find(dme.tm);           %% indices of non-zero tap ratios
            tm(i) = dme.tm(i);          %% assign non-zero tap ratios

            b = 1 ./ dme.x;             %% series susceptance
            b = b ./ tm;
            Pfinj = b .* (-dme.ta);
            obj.B = sparse( ...
                [1:nl 1:nl nl+1:2*nl nl+1:2*nl]', ...
                [1:nl nl+1:2*nl 1:nl nl+1:2*nl]', ...
                [b; -b; -b; b], ...
                2*nl, 2*nl );
            obj.p = [Pfinj + dme.g_fr; -Pfinj + dme.g_to];
        end

        %%-----  PF methods  -----
        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            %% branch active power flow
            pp = nm.get_idx('port');
            pl_fr = nm.soln.gp(pp.i1.branch(1):pp.iN.branch(1)) * dm.base_mva;
            pl_to = nm.soln.gp(pp.i1.branch(2):pp.iN.branch(2)) * dm.base_mva;

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pl_fr(dme.on) = pl_fr;
            dme.tab.pl_to(dme.on) = pl_to;
        end

        %%-----  OPF methods  -----
        function obj = opf_add_constraints(obj, mm, nm, dm, mpopt)
            %% find branches with flow limits
            dme = obj.data_model_element(dm);
            ibr = find(dme.rate_a ~= 0 & dme.rate_a < 1e10);
            nl2 = length(ibr);      %% number of constrained branches
            mm.userdata.flow_constrained_branch_idx = ibr;

            if nl2
                %% limits
                flow_max = dme.rate_a(ibr); %% RATE_A

                %% branch flow constraints
                [B, K, p] = obj.get_params(ibr);
                Af = B * obj.C';
                mm.add_lin_constraint('Pf', Af, -p-flow_max, -p+flow_max, ...
                    nm.va.order);
            end

            %% branch voltage angle difference limits
            [Aang, lang, uang, iang] = ...
                dm.elements.branch.opf_branch_ang_diff_params(...
                    dm, mpopt.opf.ignore_angle_lim);
            if length(iang)
                mm.add_lin_constraint('ang', Aang, lang, uang, {'Va'});
            end
            mm.userdata.ang_diff_constrained_branch_idx = iang;
        end

        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
            %% branch active power flow
            pp = nm.get_idx('port');
            pl_fr = nm.soln.gp(pp.i1.branch(1):pp.iN.branch(1)) * dm.base_mva;
            pl_to = nm.soln.gp(pp.i1.branch(2):pp.iN.branch(2)) * dm.base_mva;

            %% shadow prices on branch flow constraints
            ibr = mm.userdata.flow_constrained_branch_idx;
            mu_flow_fr_ub = zeros(obj.nk, 1);
            mu_flow_to_ub = mu_flow_fr_ub;
            if length(ibr)
                ll = mm.get_idx('lin');
                lambda = mm.soln.lambda;
                mu_flow_fr_ub(ibr) = lambda.mu_u(ll.i1.Pf:ll.iN.Pf);
                mu_flow_to_ub(ibr) = lambda.mu_l(ll.i1.Pf:ll.iN.Pf);
            end

            %% shadow prices on angle difference limits
            [mu_vad_lb, mu_vad_ub] = obj.opf_branch_ang_diff_prices(mm);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pl_fr(dme.on) = pl_fr;
            dme.tab.pl_to(dme.on) = pl_to;
            dme.tab.mu_flow_fr_ub(dme.on) = mu_flow_fr_ub / dm.base_mva;
            dme.tab.mu_flow_to_ub(dme.on) = mu_flow_to_ub / dm.base_mva;
            dme.tab.mu_vad_lb(dme.on) = mu_vad_lb * pi/180;
            dme.tab.mu_vad_ub(dme.on) = mu_vad_ub * pi/180;
        end
    end     %% methods
end         %% classdef
