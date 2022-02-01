classdef mp_network_acps < mp_network_acp & mp_form_acps

%   MATPOWER
%   Copyright (c) 2019-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         va = [];
%         vm = [];
%     end
    
    methods
        %%-----  CPF methods  -----
        function varargout = cpf_expand_z_warmstart(obj, ad, varargin)
            %% expand input tangent z vectors to all nodes + lambda
            varargout = cell(size(varargin));
            i = [ad.pv; ad.pq; obj.nv/2 + ad.pq; obj.nv+1];
            for k = 1:length(varargin)
                z = zeros(obj.nv, 1);
                z(i) = varargin{k};
                varargout{k} = z;
            end
        end

        function opt = cpf_solve_opts_warmstart(obj, opt, ws, mm)
            ad = mm.aux_data;

            %% update warm start states and tangent vectors
            ws.x  = [angle(ws.cV([ad.pv; ad.pq])); abs(ws.cV(ad.pq)); ws.clam];
            ws.xp = [angle(ws.pV([ad.pv; ad.pq])); abs(ws.pV(ad.pq)); ws.plam];
            opt.x0 = ws.x;   %% ignored, overridden by ws.x

            %% reduce tangent vectors for this mm
            i = [ad.pv; ad.pq; obj.nv/2 + ad.pq; obj.nv+1];
            ws.z  = ws.z(i);
            ws.zp = ws.zp(i);
            opt.warmstart = ws;
        end

        %%-----  OPF methods  -----
        function [lam_p, lam_q] = opf_node_power_balance_prices(obj, mm)
            %% shadow prices on node power balance
            nne = mm.get_idx('nle');
            lambda = mm.soln.lambda;
            lam_p = lambda.eqnonlin(nne.i1.Pmis:nne.iN.Pmis);
            lam_q = lambda.eqnonlin(nne.i1.Qmis:nne.iN.Qmis);
        end
    end     %% methods
end         %% classdef
