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

            tap = ones(nl, 1);          %% default tap ratio = 1
            i = find(dme.tap);          %% indices of non-zero tap ratios
            tap(i) = dme.tap(i);        %% assign non-zero tap ratios

            b = 1 ./ dme.X;             %% series susceptance
            b = b ./ tap;
            Pfinj = b .* (-dme.shift);
            obj.B = sparse( ...
                [1:nl 1:nl nl+1:2*nl nl+1:2*nl]', ...
                [1:nl nl+1:2*nl 1:nl nl+1:2*nl]', ...
                [b; -b; -b; b], ...
                2*nl, 2*nl );
            obj.p = [Pfinj; -Pfinj];
        end

        function add_opf_constraints(obj, nm, mm, dm, mpopt)
            %% find branches with flow limits
            dme = obj.data_model_element(dm);
            il = find(dme.rate_a ~= 0 & dme.rate_a < 1e10);
            nl2 = length(il);       %% number of constrained lines

            if nl2
                %% limits
                flow_max = dme.rate_a(il);  %% RATE_A

                %% branch flow constraints
                [B, K, p] = obj.get_params(il);
                Af = B * obj.C';
                mm.add_lin_constraint('Pf', Af, -p-flow_max, -p+flow_max, ...
                    {nm.va.order(:).name});
            end

            %% branch voltage angle difference limits
            [Aang, lang, uang, iang] = ...
                dm.branch_angle_diff_constraint(mpopt.opf.ignore_angle_lim);
            mm.add_lin_constraint('ang', Aang, lang, uang, {'Va'});
            mm.userdata.iang = iang;
        end
    end     %% methods
end         %% classdef
