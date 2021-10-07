classdef mm_pf_shared_acpi < mm_pf_shared_ac_i

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

            %% voltage angles
            st = nm.(vvars{1});
            d = st.data;
            mmx_i1 = obj.var.N + 1;
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
            if ad.npv
                obj.aux_data.var_map{end+1} = ...
                    {vvars{1}, [], [], ad.pv, mmx_i1, mmx_iN, []};
            end

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
            mmx_iN = obj.var.N;
            if ad.npq
                obj.aux_data.var_map{end+1} = ...
                    {vvars{1}, [], [], ad.pq, mmx_i1, mmx_iN, []};
            end

            %% reactive injections
            v_ = ad.vm .* exp(1j * ad.va);
            z_ = ad.zr + 1j * ad.zi;
            Qpv = nm.C(ad.pv, :) * imag( nm.port_inj_power([v_; z_], 1) );
            Qg_pv = Qpv + ad.zi(ad.zi_idx);
            mmx_i1 = obj.var.N + 1;
            obj.add_var('Qg_pv', ad.npv, Qg_pv);
            mmx_iN = obj.var.N;
            if ad.npv
                obj.aux_data.var_map{end+1} = ...
                    {'zi', [], [], ad.zi_idx, mmx_i1, mmx_iN, []};
            end

            %% voltage magnitudes
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
            mmx_iN = obj.var.N;
            if ad.npq
                obj.aux_data.var_map{end+1} = ...
                    {vvars{2}, [], [], ad.pq, mmx_i1, mmx_iN, []};
            end
        end
    end     %% methods
end         %% classdef
