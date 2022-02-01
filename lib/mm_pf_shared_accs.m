classdef mm_pf_shared_accs < mm_pf_shared_acc

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end
    
    methods
        function obj = add_pf_system_vars(obj, nm, dm, mpopt)
            %% get model variables
            vvars = nm.model_vvars();

            %% index vectors
            ad = obj.aux_data;

            %% voltage real part
            st = nm.(vvars{1});
            d = st.data;
            mmx_i1 = obj.var.N + 1;
            for k = 1:st.NS
                name = st.order(k).name;
                idx = st.order(k).idx;
                pq = ad.node_type_by_elm(k).pq;
                npq = length(pq);
                if isempty(idx)
                    obj.add_var([name '_pq'], npq, d.v0.(name)(pq), d.vl.(name)(pq), d.vu.(name)(pq));
                else
                    if all(cell2mat(idx) == 1)
                        dim = size(st.idx.N.(name));
                        if dim(end) == 1, dim(end) = []; end    %% delete trailing 1
                        obj.init_indexed_name('var', [name '_pq'], num2cell(dim));
                    end
                    sc = struct('type', {'{}', '()'}, 'subs', {idx, {pq}});
                    v0 = subsref(d.v0.(name), sc);
                    vl = subsref(d.vl.(name), sc);
                    vu = subsref(d.vu.(name), sc);
                    obj.add_var([name '_pq'], idx, npq, v0, vl, vu);
                end
            end

            for k = 1:st.NS
                name = st.order(k).name;
                idx = st.order(k).idx;
                pv = ad.node_type_by_elm(k).pv;
                npv = length(pv);
                if isempty(idx)
                    obj.add_var([name '_pv'], npv, d.v0.(name)(pv), d.vl.(name)(pv), d.vu.(name)(pv));
                else
                    if all(cell2mat(idx) == 1)
                        dim = size(st.idx.N.(name));
                        if dim(end) == 1, dim(end) = []; end    %% delete trailing 1
                        obj.init_indexed_name('var', [name '_pv'], num2cell(dim));
                    end
                    sc = struct('type', {'{}', '()'}, 'subs', {idx, {pv}});
                    v0 = subsref(d.v0.(name), sc);
                    vl = subsref(d.vl.(name), sc);
                    vu = subsref(d.vu.(name), sc);
                    obj.add_var([name '_pv'], idx, npv, v0, vl, vu);
                end
            end
            mmx_iN = obj.var.N;
            if ad.npv || ad.npq
                obj.aux_data.var_map{end+1} = ...
                    {vvars{1}, [], [], [ad.pq; ad.pv], mmx_i1, mmx_iN, []};
            end

            %% voltage imaginary part
            st = nm.(vvars{2});
            d = st.data;
            mmx_i1 = obj.var.N + 1;
            for k = 1:st.NS
                name = st.order(k).name;
                idx = st.order(k).idx;
                pq = ad.node_type_by_elm(k).pq;
                npq = length(pq);
                if isempty(idx)
                    obj.add_var([name '_pq'], npq, d.v0.(name)(pq), d.vl.(name)(pq), d.vu.(name)(pq));
                else
                    if all(cell2mat(idx) == 1)
                        dim = size(st.idx.N.(name));
                        if dim(end) == 1, dim(end) = []; end    %% delete trailing 1
                        obj.init_indexed_name('var', [name '_pq'], num2cell(dim));
                    end
                    sc = struct('type', {'{}', '()'}, 'subs', {idx, {pq}});
                    v0 = subsref(d.v0.(name), sc);
                    vl = subsref(d.vl.(name), sc);
                    vu = subsref(d.vu.(name), sc);
                    obj.add_var([name '_pq'], idx, npq, v0, vl, vu);
                end
            end

            for k = 1:st.NS
                name = st.order(k).name;
                idx = st.order(k).idx;
                pv = ad.node_type_by_elm(k).pv;
                npv = length(pv);
                if isempty(idx)
                    obj.add_var([name '_pv'], npv, d.v0.(name)(pv), d.vl.(name)(pv), d.vu.(name)(pv));
                else
                    if all(cell2mat(idx) == 1)
                        dim = size(st.idx.N.(name));
                        if dim(end) == 1, dim(end) = []; end    %% delete trailing 1
                        obj.init_indexed_name('var', [name '_pv'], num2cell(dim));
                    end
                    sc = struct('type', {'{}', '()'}, 'subs', {idx, {pv}});
                    v0 = subsref(d.v0.(name), sc);
                    vl = subsref(d.vl.(name), sc);
                    vu = subsref(d.vu.(name), sc);
                    obj.add_var([name '_pv'], idx, npv, v0, vl, vu);
                end
            end
            mmx_iN = obj.var.N;
            if ad.npv || ad.npq
                obj.aux_data.var_map{end+1} = ...
                    {vvars{2}, [], [], [ad.pq; ad.pv], mmx_i1, mmx_iN, []};
            end
        end

        function [f, J] = pf_node_balance_equations(obj, x, nm)
            %% index vector
            ad = obj.aux_data;
            pqv = [ad.pq; ad.pv];

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.pf_convert_x(x, nm, 1);

            %% incidence matrix
            C = nm.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                var_names = cellfun(@(x)x{1}, ad.var_map, 'UniformOutput', false);
                dz = any(strcmp(var_names, 'zr')) || ...
                     any(strcmp(var_names, 'zi'));
                if dz
                    [S, dS.vr, dS.vi, dS.zr, dS.zi] = nm.port_inj_power([v_; z_], 1);
                else
                    [S, dS.vr, dS.vi] = nm.port_inj_power([v_; z_], 1);
                end
                dS.vr = C * dS.vr;
                dS.vi = C * dS.vi;
                if dz
                    dS.zr = C * dS.zr;
                    dS.zi = C * dS.zi;
                end

                %% derivatives of voltage magnitudes (for PV buses)
                nn = nm.node.N;
                dV2.vr = sparse(ad.pv, ad.pv, 2*real(v_(ad.pv)), nn, nn);
                dV2.vi = sparse(ad.pv, ad.pv, 2*imag(v_(ad.pv)), nn, nn);
                dV2.zr = sparse(nn, nn);
                dV2.zi = dV2.zr;
                JJ = cell(3, length(ad.var_map));

                for k = 1:length(ad.var_map)
                    m = ad.var_map{k};
                    name = m{1};
                    if ~isempty(m{2})       %% i1:iN
                        i1 = m{2};
                        iN = m{3};
                        JJ{1, k} = real(dS.(name)(pqv,   i1:iN));
                        JJ{2, k} = imag(dS.(name)(ad.pq, i1:iN));
                        JJ{3, k} = dV2.(name)(ad.pv, i1:iN);
                    elseif isempty(m{4})    %% :
                        JJ{1, k} = real(dS.(name)(pqv,   :));
                        JJ{2, k} = imag(dS.(name)(ad.pq, :));
                        JJ{3, k} = dV2.(name)(ad.pv, :);
                    else                    %% idx
                        idx = m{4};
                        JJ{1, k} = real(dS.(name)(pqv,   idx));
                        JJ{2, k} = imag(dS.(name)(ad.pq, idx));
                        JJ{3, k} = dV2.(name)(ad.pv, idx);
                    end
                end
                J = vertcat( horzcat(JJ{1, :}), ...
                             horzcat(JJ{2, :}), ...
                             horzcat(JJ{3, :})  );
            else
                %% get port power injections (w/o derivatives)
                S = nm.port_inj_power([v_; z_], 1);
            end

            %% nodal power balance + voltage magnitude mismatch
            SS = C * S;
            vmm = v_(ad.pv) .* conj(v_(ad.pv)) - ad.vr(ad.pv).^2 - ad.vi(ad.pv).^2;
            f = [real(SS(pqv)); imag(SS(ad.pq)); vmm];
        end
    end     %% methods
end         %% classdef
