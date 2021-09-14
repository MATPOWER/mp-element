classdef mp_network_acpi < mp_network_acp & mp_form_acpi

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
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
        %%-----  PF methods  -----
        function ad = pf_aux_data(obj, dm, mpopt)
            %% call parent method
            ad = pf_aux_data@mp_network_ac(obj, dm, mpopt);

            %% build additional aux data
            N = obj.C(ad.pv, :) * obj.N;    %% z coefficients for z @ PV nodes
            [ii, jj, ~] = find(N == -1);    %% find element representing 1st
            [~, ia] = unique(ii, 'first');  %% direct injection at each PV node
            [~, ib] = sort(ii(ia));         %% sort by PV node

            %% save additional aux data
            ad.zi_idx = jj(ia(ib));         %% z-var index corresp to PV nodes
        end

        function obj = pf_add_system_vars(obj, mm, nm, dm, mpopt)
            %% get model variables
            vvars = obj.model_vvars();

            %% index vectors
            ad = mm.aux_data;

            %% voltage angles
            st = obj.(vvars{1});
            d = st.data;
            mmx_i1 = mm.var.N + 1;
            for k = 1:st.NS
                name = st.order(k).name;
                idx = st.order(k).idx;
                pv = ad.node_type_by_elm(k).pv;
                npv = length(pv);
                if isempty(idx)
                    mm.add_var([name '_pv'], npv, d.v0.(name)(pv), d.vl.(name)(pv), d.vu.(name)(pv));
                else
                    if all(cell2mat(idx) == 1)
                        dim = size(st.idx.N.(name));
                        if dim(end) == 1, dim(end) = []; end   %% delete trailing 1
                        mm.init_indexed_name('var', [name '_pv'], num2cell(dim));
                    end
                    sc = struct('type', {'{}', '()'}, 'subs', {idx, {pv}});
                    v0 = subsref(d.v0.(name), sc);
                    vl = subsref(d.vl.(name), sc);
                    vu = subsref(d.vu.(name), sc);
                    mm.add_var([name '_pv'], idx, npv, v0, vl, vu);
                end
            end
            mmx_iN = mm.var.N;
            mm.aux_data.var_map{end+1} = ...
                {vvars{1}, [], [], ad.pv, mmx_i1, mmx_iN, []};

            mmx_i1 = mm.var.N + 1;
            for k = 1:st.NS
                name = st.order(k).name;
                idx = st.order(k).idx;
                pq = ad.node_type_by_elm(k).pq;
                npq = length(pq);
                if isempty(idx)
                    mm.add_var([name '_pq'], npq, d.v0.(name)(pq), d.vl.(name)(pq), d.vu.(name)(pq));
                else
                    if all(cell2mat(idx) == 1)
                        dim = size(st.idx.N.(name));
                        if dim(end) == 1, dim(end) = []; end   %% delete trailing 1
                        mm.init_indexed_name('var', [name '_pq'], num2cell(dim));
                    end
                    sc = struct('type', {'{}', '()'}, 'subs', {idx, {pq}});
                    v0 = subsref(d.v0.(name), sc);
                    vl = subsref(d.vl.(name), sc);
                    vu = subsref(d.vu.(name), sc);
                    mm.add_var([name '_pq'], idx, npq, v0, vl, vu);
                end
            end
            mmx_iN = mm.var.N;
            mm.aux_data.var_map{end+1} = ...
                {vvars{1}, [], [], ad.pq, mmx_i1, mmx_iN, []};

            %% reactive injections
            v_ = ad.vm .* exp(1j * ad.va);
            z_ = ad.zr + 1j * ad.zi;
            Qpv = obj.C(ad.pv, :) * imag( obj.port_inj_power([v_; z_], 1) );
            Qg_pv = Qpv + ad.zi(ad.zi_idx);
            mmx_i1 = mm.var.N + 1;
            mm.add_var('Qg_pv', ad.npv, Qg_pv);
            mmx_iN = mm.var.N;
            mm.aux_data.var_map{end+1} = ...
                {'zi', [], [], ad.zi_idx, mmx_i1, mmx_iN, []};

            %% voltage magnitudes
            st = obj.(vvars{2});
            d = st.data;
            mmx_i1 = mm.var.N + 1;
            for k = 1:st.NS
                name = st.order(k).name;
                idx = st.order(k).idx;
                pq = ad.node_type_by_elm(k).pq;
                npq = length(pq);
                if isempty(idx)
                    mm.add_var([name '_pq'], npq, d.v0.(name)(pq), d.vl.(name)(pq), d.vu.(name)(pq));
                else
                    if all(cell2mat(idx) == 1)
                        dim = size(st.idx.N.(name));
                        if dim(end) == 1, dim(end) = []; end   %% delete trailing 1
                        mm.init_indexed_name('var', [name '_pq'], num2cell(dim));
                    end
                    sc = struct('type', {'{}', '()'}, 'subs', {idx, {pq}});
                    v0 = subsref(d.v0.(name), sc);
                    vl = subsref(d.vl.(name), sc);
                    vu = subsref(d.vu.(name), sc);
                    mm.add_var([name '_pq'], idx, npq, v0, vl, vu);
                end
            end
            mmx_iN = mm.var.N;
            mm.aux_data.var_map{end+1} = ...
                {vvars{2}, [], [], ad.pq, mmx_i1, mmx_iN, []};
        end

        function [f, J] = pf_node_balance_equations(obj, x, ad)
            %% index vector
            pvq = [ad.pv; ad.pq];

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.pf_convert_x(x, ad, 1);

            %% incidence matrix
            C = obj.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                [I, Iva, Ivm, Izr, Izi] = obj.port_inj_current([v_; z_], 1);

                IIva = C * Iva;
                IIvm = C * Ivm;
                IIzi = C * Izi;
                IIvm(:, ad.pv) = IIzi(:, ad.zi_idx);    %% dImis_dQg

                J = [   real(IIva(pvq, pvq))    real(IIvm(pvq, pvq));
                        imag(IIva(pvq, pvq))    imag(IIvm(pvq, pvq))    ];
            else
                %% get port power injections (w/o derivatives)
                I = obj.port_inj_current([v_; z_], 1);
            end

            %% nodal power balance
            II = C * I;
            f = [real(II(pvq)); imag(II(pvq))];
        end

        function obj = pf_add_node_balance_constraints(obj, mm, dm, mpopt)
            %% power balance constraints
            ad = mm.aux_data;
            npvq = ad.npv+ad.npq;
            fcn = @(x)pf_node_balance_equations(obj, x, ad);
            mm.add_nln_constraint({'Irmis', 'Iimis'}, [npvq;npvq], 1, fcn, []);
        end


        %%-----  OPF methods  -----
        function opf_add_node_balance_constraints(obj, mm)
            %% power balance constraints
            nn = obj.node.N;            %% number of nodes
            fcn_mis = @(x)opf_current_balance_fcn(obj, obj.opf_convert_x(x, mm.aux_data));
            hess_mis = @(x, lam)opf_current_balance_hess(obj, ...
                obj.opf_convert_x(x, mm.aux_data), lam);
            mm.add_nln_constraint({'rImis', 'iImis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end

        function [lam_p, lam_q] = opf_node_power_balance_prices(obj, mm)
            %% shadow prices on node power balance
            nne = mm.get_idx('nle');

            %% convert current balance shadow prices to equivalent lam_p and lam_q
            %% P + jQ = (Vr + jVi) * (M - jN)
            %% M = (Vr P + Vi Q) / (Vr^2 + Vi^2)
            %% N = (Vi P - Vr Q) / (Vr^2 + Vi^2)
            %% lam_p = df/dP = df/dM * dM/dP + df/dN + dN/dP
            %% lam_q = df/dQ = df/dM * dM/dQ + df/dN + dN/dQ
            V = obj.soln.v;
            lambda = mm.soln.lambda;
            VV = V ./ (V .* conj(V));   %% V / vm^2
            VVr = real(VV);
            VVi = imag(VV);
            lamM = lambda.eqnonlin(nne.i1.rImis:nne.iN.rImis);
            lamN = lambda.eqnonlin(nne.i1.iImis:nne.iN.iImis);
            lam_p = (VVr.*lamM + VVi.*lamN);
            lam_q = (VVi.*lamM - VVr.*lamN);
        end
    end     %% methods
end         %% classdef
