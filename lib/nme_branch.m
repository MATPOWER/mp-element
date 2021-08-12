classdef nme_branch < nm_element

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
        %% constructor
        function obj = nme_branch()
            obj@nm_element();
            obj.name = 'branch';
            obj.np = 2;             %% this is a 2 port element
        end

        function [fidx, tidx] = node_indices(obj, nm, dm, dme)
            bus_dme = dm.elements.bus;
            nidx = nm.get_node_idx('bus');  %% node indices for 'bus'
            f = dme.fbus(dme.on);
            t = dme.tbus(dme.on);
            fbidx = bus_dme.i2on(f);        %% online bus indices branch "from"
            tbidx = bus_dme.i2on(t);        %% online bus indices branch "to"
            fidx = nidx(fbidx);             %% branch "from" port node indices
            tidx = nidx(tbidx);             %% branch "to" port node indices
        end

        function obj = build_params(obj, nm, dm)
            dme = obj.data_model_element(dm);

            %% incidence matrices
            [fidx, tidx] = obj.node_indices(nm, dm, dme);
            obj.C = obj.incidence_matrix(nm.getN('node'), fidx, tidx);
            obj.D = obj.incidence_matrix(nm.getN('state'));
        end

        function [mu_vad_lb, mu_vad_ub] = opf_branch_ang_diff_prices(obj, mm)
            %% shadow prices on angle difference limits
            iang = mm.userdata.ang_diff_constrained_branch_idx;
            mu_vad_lb = zeros(obj.nk, 1);
            mu_vad_ub = mu_vad_lb;
            if length(iang)
                ll = mm.get_idx('lin');
                lambda = mm.soln.lambda;
                mu_vad_lb(iang) = lambda.mu_l(ll.i1.ang:ll.iN.ang);
                mu_vad_ub(iang) = lambda.mu_u(ll.i1.ang:ll.iN.ang);
            end
        end
    end     %% methods
end         %% classdef
