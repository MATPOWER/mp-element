classdef mme_buslink_pf_ac < mme_buslink

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'buslink';
%     end
    
    methods
        function obj = add_vars(obj, mm, nm, dm, mpopt)
            nme = obj.network_model_element(nm);

            mm.init_indexed_name('var', 'Plink', {3});
            mm.init_indexed_name('var', 'Qlink', {3});
            mmx_i1 = mm.var.N + 1;
            for p = 1:3
                mm.add_var('Plink', {p}, nme.nk, 0, -Inf, Inf);
            end
            mmx_iN = mm.var.N;
            mm.aux_data.var_map{end+1} = ...
                {'zr', nm.zr.idx.i1.Plink(1), nm.zr.idx.iN.Plink(end), [], mmx_i1, mmx_iN, []};

            mmx_i1 = mm.var.N + 1;
            for p = 1:3
                mm.add_var('Qlink', {p}, nme.nk, 0, -Inf, Inf);
            end
            mmx_iN = mm.var.N;
            mm.aux_data.var_map{end+1} = ...
                {'zi', nm.zi.idx.i1.Qlink(1), nm.zi.idx.iN.Qlink(end), [], mmx_i1, mmx_iN, []};
        end

        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            nme = obj.network_model_element(nm);

            %% add Qlink constraints for PV node on 3-phase side

            %% indices of buslinks connected to PV nodes via port 2 (3&4)
            ad = mm.aux_data;
            [~, k, ~] = find(nme.C(ad.pv, nme.nk+1:2*nme.nk));
            n = length(k);
            if n
                I = sparse(1:n, k, 1, n, nme.nk);
                zz = sparse(n, nme.nk);
                A = [I -I zz; I zz -I];
                b = zeros(2*n, 1);
                vs = struct('name', {'Qlink', 'Qlink', 'Qlink'}, ...
                            'idx', {{1}, {2}, {3}} );
                mm.add_lin_constraint('buslink_qlink', A, b, b, vs);
            end
        end
    end     %% methods
end         %% classdef
