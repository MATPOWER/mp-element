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

        function add_opf_constraints(obj, nm, mm, dm, mpopt)
            %% find branches with flow limits
            dme = obj.data_model_element(dm);
            il = find(dme.rate_a ~= 0 & dme.rate_a < 1e10);
            nl2 = length(il);       %% number of constrained lines

            if nl2
                %% port indexes
                nl = obj.nk;
                idx = [il; nl+il];

                %% limits
                flow_max = dme.rate_a(il);  %% RATE_A

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
                    error('nme_branch_ac/add_opf_constraints: MPOPT.opf.flow_lim = ''%s'' not yet implemented.', mpopt.opf.flow_lim);
                end
            
                mm.add_nln_constraint({'Sf', 'St'}, [nl2;nl2], 0, fcn_flow, hess_flow);
            end
        end
    end     %% methods
end         %% classdef
