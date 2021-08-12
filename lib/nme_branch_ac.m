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

            tap = ones(nl, 1);          %% default tap ratio = 1
            i = find(dme.tap);          %% indices of non-zero tap ratios
            tap(i) = dme.tap(i);        %% assign non-zero tap ratios
            tap = tap .* exp(1j * dme.shift);   %% add phase shifters

            Ys = 1 ./ (dme.R + 1j * dme.X);     %% series admittance
            Bc = dme.B;                         %% line charging susceptance
            Ytt = Ys + 1j*Bc/2;
            Yff = Ytt ./ (tap .* conj(tap));
            Yft = - Ys ./ conj(tap);
            Ytf = - Ys ./ tap;

            obj.Y = sparse( ...
                [1:nl 1:nl nl+1:2*nl nl+1:2*nl]', ...
                [1:nl nl+1:2*nl 1:nl nl+1:2*nl]', ...
                [Yff; Yft; Ytf; Ytt], 2*nl, 2*nl );
        end

        %%-----  PF methods  -----
        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            %% branch active power flow
            pp = nm.get_idx('port');
            Sf = nm.soln.gs_(pp.i1.branch(1):pp.iN.branch(1)) * dm.baseMVA;
            St = nm.soln.gs_(pp.i1.branch(2):pp.iN.branch(2)) * dm.baseMVA;

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pl_fr(dme.on) = real(Sf);
            dme.tab.ql_fr(dme.on) = imag(Sf);
            dme.tab.pl_to(dme.on) = real(St);
            dme.tab.ql_to(dme.on) = imag(St);
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
                        nm.opf_convert_x(x), nm, idx, ...
                        [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_apparent_power_lim_hess(obj, ...
                        nm.opf_convert_x(x), lam, nm, idx);
                elseif lim_type == 'P'
                    fcn_flow = @(x)port_active_power_lim_fcn(obj, ...
                        nm.opf_convert_x(x), nm, idx, [flow_max; flow_max]);
                    hess_flow = @(x, lam)port_active_power_lim_hess(obj, ...
                        nm.opf_convert_x(x), lam, nm, idx);
                elseif lim_type == '2' || lim_type == 'P'
                    fcn_flow = @(x)port_active_power2_lim_fcn(obj, ...
                        nm.opf_convert_x(x), nm, idx, ...
                        [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_active_power2_lim_hess(obj, ...
                        nm.opf_convert_x(x), lam, nm, idx);
                elseif lim_type == 'I'
                    fcn_flow = @(x)port_current_lim_fcn(obj, ...
                        nm.opf_convert_x(x), nm, idx, ...
                        [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_current_lim_hess(obj, ...
                        nm.opf_convert_x(x), lam, nm, idx);
                else
                    error('nme_branch_ac/opf_add_constraints: MPOPT.opf.flow_lim = ''%s'' not yet implemented.', mpopt.opf.flow_lim);
                end
            
                mm.add_nln_constraint({'Sf', 'St'}, [nl2;nl2], 0, fcn_flow, hess_flow);
            end
        end

        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
            %% branch active power flow
            pp = nm.get_idx('port');
            Sf = nm.soln.gs_(pp.i1.branch(1):pp.iN.branch(1)) * dm.baseMVA;
            St = nm.soln.gs_(pp.i1.branch(2):pp.iN.branch(2)) * dm.baseMVA;

            %% shadow prices on branch flow constraints
            ibr = mm.userdata.flow_constrained_branch_idx;
            muSf = zeros(obj.nk, 1);
            muSt = muSf;
            if length(ibr)
                lim_type = upper(mpopt.opf.flow_lim(1));
                nni = mm.get_idx('nli');
                lambda = mm.soln.lambda;
                if lim_type == 'P'
                    muSf(ibr) = lambda.ineqnonlin(nni.i1.Sf:nni.iN.Sf);
                    muSt(ibr) = lambda.ineqnonlin(nni.i1.St:nni.iN.St);
                else
                    rate_a = obj.data_model_element(dm).rate_a(ibr);
                    muSf(ibr) = 2 * lambda.ineqnonlin(nni.i1.Sf:nni.iN.Sf) .* rate_a;
                    muSt(ibr) = 2 * lambda.ineqnonlin(nni.i1.St:nni.iN.St) .* rate_a;
                end
            end

            %% shadow prices on angle difference limits
            [muAngmin, muAngmax] = obj.opf_branch_ang_diff_prices(mm);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.tab.pl_fr(dme.on) = real(Sf);
            dme.tab.ql_fr(dme.on) = imag(Sf);
            dme.tab.pl_to(dme.on) = real(St);
            dme.tab.ql_to(dme.on) = imag(St);
            dme.tab.mu_flow_fr_ub(dme.on) = muSf / dm.baseMVA;
            dme.tab.mu_flow_to_ub(dme.on) = muSt / dm.baseMVA;
            dme.tab.mu_vad_lb(dme.on) = muAngmin * pi/180;
            dme.tab.mu_vad_ub(dme.on) = muAngmax * pi/180;
        end
    end     %% methods
end         %% classdef
