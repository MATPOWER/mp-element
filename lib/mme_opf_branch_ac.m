classdef mme_opf_branch_ac < mme_branch

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'branch';
%     end
    
    methods
        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            %% find branches with flow limits
            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);
            ibr = find(dme.rate_a ~= 0 & dme.rate_a < 1e10);
            nl2 = length(ibr);      %% number of constrained branches
            mm.userdata.flow_constrained_branch_idx = ibr;

            if nl2
                %% port indexes
                nl = nme.nk;
                idx = [ibr; nl+ibr];

                %% limits
                flow_max = dme.rate_a(ibr); %% RATE_A

                %% branch flow constraints
                lim_type = upper(mpopt.opf.flow_lim(1));
                if lim_type == 'S'
                    fcn_flow = @(x)port_apparent_power_lim_fcn(nme, ...
                        nm.opf_convert_x(x, mm.aux_data), nm, idx, ...
                        [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_apparent_power_lim_hess(nme, ...
                        nm.opf_convert_x(x, mm.aux_data), lam, nm, idx);
                elseif lim_type == 'P'
                    fcn_flow = @(x)port_active_power_lim_fcn(nme, ...
                        nm.opf_convert_x(x, mm.aux_data), nm, idx, [flow_max; flow_max]);
                    hess_flow = @(x, lam)port_active_power_lim_hess(nme, ...
                        nm.opf_convert_x(x, mm.aux_data), lam, nm, idx);
                elseif lim_type == '2' || lim_type == 'P'
                    fcn_flow = @(x)port_active_power2_lim_fcn(nme, ...
                        nm.opf_convert_x(x, mm.aux_data), nm, idx, ...
                        [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_active_power2_lim_hess(nme, ...
                        nm.opf_convert_x(x, mm.aux_data), lam, nm, idx);
                elseif lim_type == 'I'
                    fcn_flow = @(x)port_current_lim_fcn(nme, ...
                        nm.opf_convert_x(x, mm.aux_data), nm, idx, ...
                        [flow_max; flow_max] .^ 2);
                    hess_flow = @(x, lam)port_current_lim_hess(nme, ...
                        nm.opf_convert_x(x, mm.aux_data), lam, nm, idx);
                else
                    error('mme_opf_branch_ac/add_constraints: MPOPT.opf.flow_lim = ''%s'' not yet implemented.', mpopt.opf.flow_lim);
                end

                mm.add_nln_constraint({'Sf', 'St'}, [nl2;nl2], 0, fcn_flow, hess_flow);
            end
        end
    end     %% methods
end         %% classdef
