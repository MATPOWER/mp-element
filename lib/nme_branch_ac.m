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
    end     %% methods
end         %% classdef
