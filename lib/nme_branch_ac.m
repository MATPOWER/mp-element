classdef nme_branch_ac < nme_branch% & mp_form_ac

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
            tm(i) = dme.tm(i);              %% assign non-zero tap ratios
            T = tm .* exp(1j * dme.ta);     %% add phase shifters

            ys = 1 ./ (dme.r + 1j * dme.x);     %% series admittance
            yff = (ys + dme.g_fr + 1j * dme.b_fr) ./ (T .* conj(T));
            ytt = ys + dme.g_to + 1j * dme.b_to;
            yft = - ys ./ conj(T);
            ytf = - ys ./ T;

            obj.Y = sparse( ...
                [1:nl 1:nl nl+1:2*nl nl+1:2*nl]', ...
                [1:nl nl+1:2*nl 1:nl nl+1:2*nl]', ...
                [yff; yft; ytf; ytt], 2*nl, 2*nl );

% same as:
%             Yff = spdiags(yff, 0, nl, nl);
%             Yft = spdiags(yft, 0, nl, nl);
%             Ytf = spdiags(ytf, 0, nl, nl);
%             Ytt = spdiags(ytt, 0, nl, nl);
%             obj.Y = [ Yff Yft;
%                       Ytf Ytt  ];
        end

        %%-----  PF methods  -----
        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            %% branch active power flow
            pp = nm.get_idx('port');
            S_fr = nm.soln.gs_(pp.i1.branch(1):pp.iN.branch(1)) * dm.base_mva;
            S_to = nm.soln.gs_(pp.i1.branch(2):pp.iN.branch(2)) * dm.base_mva;

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pl_fr(dme.on) = real(S_fr);
            dme.tab.ql_fr(dme.on) = imag(S_fr);
            dme.tab.pl_to(dme.on) = real(S_to);
            dme.tab.ql_to(dme.on) = imag(S_to);
        end

        %%-----  OPF methods  -----
        function obj = opf_add_constraints(obj, mm, nm, dm, mpopt)
            %% find branches with flow limits
            dme = obj.data_model_element(dm);
            ibr = find(dme.rate_a ~= 0 & dme.rate_a < 1e10);
            nl2 = length(ibr);      %% number of constrained branches
            mm.userdata.flow_constrained_branch_idx = ibr;

            if nl2
                %% port indexes
                nl = obj.nk;
                idx = [ibr; nl+ibr];

                %% limits
                flow_max = dme.rate_a(ibr); %% RATE_A

                %% branch flow constraints
                lim_type = upper(mpopt.opf.flow_lim(1));
                if lim_type == 'S'
                    fcn_flow = @(x)port_apparent_power_lim_fcn(obj, ...
                        nm.opf_convert_x(x, mm.aux_data), nm, idx, ...
                        [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_apparent_power_lim_hess(obj, ...
                        nm.opf_convert_x(x, mm.aux_data), lam, nm, idx);
                elseif lim_type == 'P'
                    fcn_flow = @(x)port_active_power_lim_fcn(obj, ...
                        nm.opf_convert_x(x, mm.aux_data), nm, idx, [flow_max; flow_max]);
                    hess_flow = @(x, lam)port_active_power_lim_hess(obj, ...
                        nm.opf_convert_x(x, mm.aux_data), lam, nm, idx);
                elseif lim_type == '2' || lim_type == 'P'
                    fcn_flow = @(x)port_active_power2_lim_fcn(obj, ...
                        nm.opf_convert_x(x, mm.aux_data), nm, idx, ...
                        [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_active_power2_lim_hess(obj, ...
                        nm.opf_convert_x(x, mm.aux_data), lam, nm, idx);
                elseif lim_type == 'I'
                    fcn_flow = @(x)port_current_lim_fcn(obj, ...
                        nm.opf_convert_x(x, mm.aux_data), nm, idx, ...
                        [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_current_lim_hess(obj, ...
                        nm.opf_convert_x(x, mm.aux_data), lam, nm, idx);
                else
                    error('nme_branch_ac/opf_add_constraints: MPOPT.opf.flow_lim = ''%s'' not yet implemented.', mpopt.opf.flow_lim);
                end
            
                mm.add_nln_constraint({'Sf', 'St'}, [nl2;nl2], 0, fcn_flow, hess_flow);
            end
        end

        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
            %% branch active power flow
            pp = nm.get_idx('port');
            S_fr = nm.soln.gs_(pp.i1.branch(1):pp.iN.branch(1)) * dm.base_mva;
            S_to = nm.soln.gs_(pp.i1.branch(2):pp.iN.branch(2)) * dm.base_mva;

            %% shadow prices on branch flow constraints
            ibr = mm.userdata.flow_constrained_branch_idx;
            mu_flow_fr_ub = zeros(obj.nk, 1);
            mu_flow_to_ub = mu_flow_fr_ub;
            if length(ibr)
                lim_type = upper(mpopt.opf.flow_lim(1));
                nni = mm.get_idx('nli');
                lambda = mm.soln.lambda;
                if lim_type == 'P'
                    mu_flow_fr_ub(ibr) = lambda.ineqnonlin(nni.i1.Sf:nni.iN.Sf);
                    mu_flow_to_ub(ibr) = lambda.ineqnonlin(nni.i1.St:nni.iN.St);
                else
                    rate_a = obj.data_model_element(dm).rate_a(ibr);
                    mu_flow_fr_ub(ibr) = 2 * lambda.ineqnonlin(nni.i1.Sf:nni.iN.Sf) .* rate_a;
                    mu_flow_to_ub(ibr) = 2 * lambda.ineqnonlin(nni.i1.St:nni.iN.St) .* rate_a;
                end
            end

            %% shadow prices on angle difference limits
            [mu_vad_lb, mu_vad_ub] = obj.opf_branch_ang_diff_prices(mm);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pl_fr(dme.on) = real(S_fr);
            dme.tab.ql_fr(dme.on) = imag(S_fr);
            dme.tab.pl_to(dme.on) = real(S_to);
            dme.tab.ql_to(dme.on) = imag(S_to);
            dme.tab.mu_flow_fr_ub(dme.on) = mu_flow_fr_ub / dm.base_mva;
            dme.tab.mu_flow_to_ub(dme.on) = mu_flow_to_ub / dm.base_mva;
            dme.tab.mu_vad_lb(dme.on) = mu_vad_lb * pi/180;
            dme.tab.mu_vad_ub(dme.on) = mu_vad_ub * pi/180;
        end
    end     %% methods
end         %% classdef
