classdef mp_network_acps < mp_network_acp & mp_form_acps

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

            switch mpopt.pf.alg
                case 'GS'
                    ad.Y = obj.C * obj.get_params([], 'Y') * obj.C';
                case 'ZG'
                    pvq = [ad.pv; ad.pq];
                    Y = obj.C * obj.get_params([], 'Y') * obj.C';
                    Y21 = Y(pvq, ad.ref);
                    Y22 = Y(pvq, pvq);
                    [L, U, p, q] = lu(Y22, 'vector');
                    [junk, iq] = sort(q);
                    [ad.Y, ad.Y21, ad.L, ad.U, ad.p, ad.iq] = ...
                        deal(Y, Y21, L, U, p, iq);
            end
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
            if ad.npv
                mm.aux_data.var_map{end+1} = ...
                    {vvars{1}, [], [], ad.pv, mmx_i1, mmx_iN, []};
            end

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
            if ad.npq
                mm.aux_data.var_map{end+1} = ...
                    {vvars{1}, [], [], ad.pq, mmx_i1, mmx_iN, []};
            end

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
            if ad.npq
                mm.aux_data.var_map{end+1} = ...
                    {vvars{2}, [], [], ad.pq, mmx_i1, mmx_iN, []};
            end
        end

        function [f, J] = pf_node_balance_equations(obj, x, ad, fdpf)
            %% index vector
            pvq = [ad.pv; ad.pq];

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.pf_convert_x(x, ad, 1);

            %% incidence matrix
            C = obj.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                var_names = cellfun(@(x)x{1}, ad.var_map, 'UniformOutput', false);
                dz = any(strcmp(var_names, 'zr')) || ...
                     any(strcmp(var_names, 'zi'));
                if dz
                    [S, dS.va, dS.vm, dS.zr, dS.zi] = obj.port_inj_power([v_; z_], 1);
                else
                    [S, dS.va, dS.vm] = obj.port_inj_power([v_; z_], 1);
                end
                dS.va = C * dS.va;
                dS.vm = C * dS.vm;
                if dz
                    dS.zr = C * dS.zr;
                    dS.zi = C * dS.zi;
                end
                JJ = cell(2, length(ad.var_map));

                for k = 1:length(ad.var_map)
                    m = ad.var_map{k};
                    name = m{1};
                    if ~isempty(m{2})       %% i1:iN
                        i1 = m{2};
                        iN = m{3};
                        JJ{1, k} = real(dS.(name)(pvq,   i1:iN));
                        JJ{2, k} = imag(dS.(name)(ad.pq, i1:iN));
                    elseif isempty(m{4})    %% :
                        JJ{1, k} = real(dS.(name)(pvq,   :));
                        JJ{2, k} = imag(dS.(name)(ad.pq, :));
                    else                    %% idx
                        idx = m{4};
                        JJ{1, k} = real(dS.(name)(pvq,   idx));
                        JJ{2, k} = imag(dS.(name)(ad.pq, idx));
                    end
                end
                J = vertcat( horzcat(JJ{1, :}), ...
                             horzcat(JJ{2, :})  );
            else
                %% get port power injections (w/o derivatives)
                S = obj.port_inj_power([v_; z_], 1);
            end

            %% nodal power balance
            if nargin > 3 && fdpf
                SS = C * S ./ abs(v_);  %% for fast-decoupled formulation
            else
                SS = C * S;
            end
            f = [real(SS(pvq)); imag(SS(ad.pq))];
        end

        function obj = pf_add_node_balance_constraints(obj, mm, dm, mpopt)
            alg = mpopt.pf.alg;
            ad = mm.aux_data;
            
            %% power balance constraints
            switch alg
                case  {'FDXB', 'FDBX'}
                    fcn = @(x)pf_node_balance_equations(obj, x, ad, 1);
                otherwise
                    fcn = @(x)pf_node_balance_equations(obj, x, ad);
            end
            mm.add_nln_constraint({'Pmis', 'Qmis'}, [ad.npv+ad.npq;ad.npq], 1, fcn, []);
        end

        function JJ = pf_fd_jac_approx(obj, mm, dm, mpopt)
            alg = mpopt.pf.alg;

            %% create copies of data model for building B prime, B double prime
            [dm1, dm2] = dm.fdpf_B_matrix_models(alg);

            %% build network models and get admittance matrices
            nm1 = feval(class(obj)).build(dm1);
            nm2 = feval(class(obj)).build(dm2);
            [Y1, L, M] = nm1.get_params([], {'Y', 'L', 'M'});
            Y2 = nm2.get_params();
            if any(any(L)) || any(any(M))
                error('mp_network_acps/df_jac_approx: fast-decoupled Jacobian approximation not implemented for models with non-zero L and/or M matrices.')
            end

            %% form reduced Bp and Bpp matrices
            ad = mm.aux_data;
            Cp  = nm1.C([ad.pv; ad.pq], :);
            Cpp = nm2.C(ad.pq, :);
            Bp  = -imag( Cp  * Y1 * Cp' );
            Bpp = -imag( Cpp * Y2 * Cpp' );
            JJ = {Bp, Bpp};
        end

        function x = pf_gs_x_update(obj, x, f, mm, dm, mpopt);
            alg = mpopt.pf.alg;
            ad = mm.aux_data;

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.pf_convert_x(x, ad, 1);

            [pv, pq, npv, npq, Y] = deal(ad.pv, ad.pq, ad.npv, ad.npq, ad.Y);
            
            %% total nodal complex bus power extractions
            SS = zeros(size(v_));
            SS(pv) = f(1:npv);
            SS(pq) = f(npv+1:npv+npq) + 1j * f(npv+npq+1:npv+2*npq);
%            SS = C * obj.port_inj_power([v_; z_], 1);

            %% complex net nodal injection (from all but constant Z elements)
            S0 = v_ .* conj(Y * v_) - SS;

            %% update voltage
            %% at PQ buses
            for k = pq'
                v_(k) = v_(k) + (conj(S0(k)/v_(k)) - Y(k,:) * v_) / Y(k, k);
            end

            %% at PV buses
            if npv
                for k = pv'
                    S0(k) = real(S0(k)) + 1j * imag( v_(k) * conj(Y(k,:) * v_) );
                    v_(k) = v_(k) + (conj(S0(k)/v_(k)) - Y(k,:) * v_) / Y(k, k);
                end
            end

            x = [angle(v_([pv; pq])); abs(v_(pq))];
        end

        function x = pf_zg_x_update(obj, x, f, mm, dm, mpopt);
            alg = mpopt.pf.alg;
            ad = mm.aux_data;

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.pf_convert_x(x, ad, 1);

            [pv, pq, ref, npv, npq] = deal(ad.pv, ad.pq, ad.ref, ad.npv, ad.npq);
            pvq = [pv; pq];

            %% build and cache S0 and Y21_v1
            if isfield(ad, 'S0')
                [Y21_v1, S0] = deal(ad.Y21_v1, ad.S0);
            else
                Y21_v1 = ad.Y21 * v_(ref);
                SS = zeros(size(v_));
                SS(pv) = f(1:npv);
                SS(pq) = f(npv+1:npv+npq) + 1j * f(npv+npq+1:npv+2*npq);
    %            SS = C * obj.port_inj_power([v_; z_], 1);

                %% complex net nodal injection (from all but constant Z elements)
                S0 = v_ .* conj(ad.Y * v_) - SS;
                S0(ref) = 0;

                %% cache 'em
                [ad.Y21_v1, ad.S0] = deal(Y21_v1, S0);
                mm.aux_data = ad;
            end

            if npv  %% update Q injections at PV buses based on vm mismatch
                %% compute and cache initial vm at PV buses and
                %% factored fast-decoupled Bpp matrix
                if isfield(ad, 'vmpv0')
                    [vmpv, vmpv0, Bpp, LBpp, UBpp, pBpp, iqBpp] = ...
                        deal(ad.vmpv, ad.vmpv0, ad.Bpp, ad.LBpp, ad.UBpp, ad.pBpp, ad.iqBpp);
                else
                    vmpv0 = abs(v_(pv));
                    vmpv = vmpv0;

                    %% modify data model to form Bpp (B double prime)
                    dm2 = dm.fdpf_B_matrix_models('FDBX');

                    %% build network models and get admittance matrices
                    nm = feval(class(obj)).build(dm2);
                    [Y2, L, M] = nm.get_params([], {'Y', 'L', 'M'});
                    if any(any(L)) || any(any(M))
                        error('mp_network_acps/pf_zg_x_update: B matrix for Z-bus Gauss w/PV buses not implemented for models with non-zero L and/or M matrices.')
                    end
                    Bpp = -nm.C * imag(Y2) * nm.C';

                    [LBpp, UBpp, pBpp, qBpp] = lu(Bpp(pq, pq), 'vector');
                    [junk, iqBpp] = sort(qBpp);

                    %% cache 'em
                    [ad.vmpv, ad.vmpv0, ad.Bpp, ad.LBpp, ad.UBpp, ad.pBpp, ad.iqBpp] = ...
                        deal(vmpv, vmpv0, Bpp, LBpp, UBpp, pBpp, iqBpp);
                    mm.aux_data = ad;
                end

                %% compute voltage mismatches at PV buses
                v_(pv) = vmpv .* v_(pv) ./ abs(v_(pv));
                dV = vmpv0 - vmpv;
%                 [max_dV, k] = max(abs(dV));
%                 fprintf('       %10.3e', max_dV)

                %% compute Q injection at current V
                %% (sometimes improves convergence)
                Qpv = imag( v_(pv) .* conj(ad.Y(pv, :) * v_) );
                S0(pv) = S0(pv) + 1j * (Qpv - imag(S0(pv)));

                % dVpq = Bpp(pq, pq) \ (-Bpp(pq, pv) * dV);
                dVpq = UBpp \  (LBpp \ (-Bpp(pq(pBpp), pv) * dV));
                dVpq = dVpq(iqBpp);
                dQ = Bpp(pv, pq) * dVpq + Bpp(pv, pv) * dV;

                %% update S0
                S0(pv) = S0(pv) + 1j * dQ;
            end

            %% complex current injections
            I2 = conj(S0(pvq) ./ v_(pvq));

            V2 = ad.U \  (ad.L \ (I2(ad.p) - Y21_v1(ad.p)));
            V2 = V2(ad.iq);

            v_(pv) = V2(1:npv);
            v_(pq) = V2(npv+1:npv+npq);
            mm.aux_data.vmpv = abs(v_(pv));

            x(1:ad.npv+ad.npq) = angle(v_(pvq));
            x(ad.npv+ad.npq+1:ad.npv+2*ad.npq) = abs(v_(pq));
        end

        %%-----  CPF methods  -----
        function cpf_add_node_balance_constraints(obj, mm, dm, mpopt)
            ad = mm.aux_data;
            
            %% continuation power balance constraints
            fcn = @(x)cpf_node_balance_equations(obj, x, ad);
            mm.add_nln_constraint({'Pmis', 'Qmis'}, [ad.npv+ad.npq;ad.npq], 1, fcn, []);
        end

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
        function opf_add_node_balance_constraints(obj, mm)
            %% power balance constraints
            nn = obj.node.N;            %% number of nodes
            fcn_mis = @(x)opf_power_balance_fcn(obj, obj.opf_convert_x(x, mm.aux_data));
            hess_mis = @(x, lam)opf_power_balance_hess(obj, ...
                obj.opf_convert_x(x, mm.aux_data), lam);
            mm.add_nln_constraint({'Pmis', 'Qmis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end

        function [lam_p, lam_q] = opf_node_power_balance_prices(obj, mm)
            %% shadow prices on node power balance
            nne = mm.get_idx('nle');
            lambda = mm.soln.lambda;
            lam_p = lambda.eqnonlin(nne.i1.Pmis:nne.iN.Pmis);
            lam_q = lambda.eqnonlin(nne.i1.Qmis:nne.iN.Qmis);
        end
    end     %% methods
end         %% classdef
