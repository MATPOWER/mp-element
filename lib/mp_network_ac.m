classdef mp_network_ac < mp_network% & mp_form_ac
%MP_NETWORK_AC Abstract class, explicitly a subclass of MP_NETWORK and
%              implicitly assumed to be a subclass of MP_FORM_AC as well

%   MATPOWER
%   Copyright (c) 2019-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        zr = [];
        zi = [];
        inln_list = {};         %% private: list of indexes of nme's w/inln
        snln_list = {};         %% private: list of indexes of nme's w/snln
        inln_hess_list = {};    %% private: list of indexes of nme's w/inln_hess
        snln_hess_list = {};    %% private: list of indexes of nme's w/snln_hess
    end

    methods
        function obj = def_set_types(obj)
            def_set_types@mp_network(obj);      %% call parent first
            obj.set_types.zr = 'NON-VOLTAGE VARS REAL (zr)';
            obj.set_types.zi = 'NON-VOLTAGE VARS IMAG (zi)';
        end

        function obj = build_params(obj, nm, dm)
            %% call parent to build individual element parameters
            build_params@mp_network(obj, nm, dm);

            %% aggregate parameters from individual elements
            obj.Y = obj.stack_matrix_params('Y', 1);
            obj.L = obj.stack_matrix_params('L', 0);
            obj.M = obj.stack_matrix_params('M', 1);
            obj.N = obj.stack_matrix_params('N', 0);
            obj.i = obj.stack_vector_params('i');
            obj.s = obj.stack_vector_params('s');

            %% add general nonlinear function if any element has one defined
            for k = 1:length(obj.elements)
                if ~isempty(obj.elements{k}.inln)
                    obj.inln_list{end+1} = k;
                    if ~isempty(obj.elements{k}.inln_hess)
                        obj.inln_hess_list{end+1} = k;
                    end
                end
                if ~isempty(obj.elements{k}.snln)
                    obj.snln_list{end+1} = k;
                    if ~isempty(obj.elements{k}.snln_hess)
                        obj.snln_hess_list{end+1} = k;
                    end
                end
            end
            if ~isempty(obj.inln_list)
                obj.inln = @(x_, sysx, idx)port_inj_nln(obj, 'i', x_, sysx, idx);
                if ~isempty(obj.inln_hess_list)
                    obj.inln_hess = @(x_, lam, sysx, idx)port_inj_nln_hess(obj, 'i', x_, lam, sysx, idx);
                end
            end
            if ~isempty(obj.snln_list)
                obj.snln = @(x_, sysx, idx)port_inj_nln(obj, 's', x_, sysx, idx);
                if ~isempty(obj.snln_hess_list)
                    obj.snln_hess = @(x_, lam, sysx, idx)port_inj_nln_hess(obj, 's', x_, lam, sysx, idx);
                end
            end
        end

        function [g, gv1, gv2, gzr, gzi] = port_inj_nln(obj, si, x_, sysx, idx)
            if nargin < 5
                idx = [];
                if nargin < 4
                    sysx = 1;
                end
            end

            %% current or power
            fcn = [si 'nln'];
            fcn_list = [fcn '_list'];

            %% initialize
            if isempty(idx)
                sel = 0;        %% all ports
                np = obj.np;
            else
                sel = 1;        %% selected ports only
                np = length(idx);
            end
            nv = obj.get_nv_(sysx);
            nz = obj.nz;
            nc = size(x_, 2);   %% num of cols in x_, for evaluating multiple x_
            g = zeros(np, nc);
            gv1 = sparse(np, nv);
            gv2 = sparse(np, nv);
            gzr = sparse(np, nz);
            gzi = sparse(np, nz);

            %% loop through elements w/gen nonlin fcns, evaluate them
            pp = obj.get_idx('port');
            if ~sysx
                ss = obj.get_idx('state');
            end
            for kk = obj.(fcn_list)
                k = kk{1};      %% index into obj.elements
                nme = obj.elements{k};
                i1 = pp.i1.(nme.name)(1);
                iN = pp.iN.(nme.name)(end);

                %% set up port index vector for nme
                if sel
                    apidx = find(idx >= i1 & idx <= iN);   %% aggregate port indices in range
                    if isempty(apidx)  %% skip if selected ports, but none in range
                        continue;
                    end
                    nme_idx = idx(apidx) - i1 + 1;  %% port index vector for nme
                else
                    nme_idx = [];                   %% all ports for nme
                end

                %% set up proper x_ for nme
                if sysx
                    nme_x_ = x_;
                else
                    if isfield(ss.i1, nme.name)
                        j1 = ss.i1.(nme.name)(1);
                        jN = ss.iN.(nme.name)(end);
                    else
                        j1 = 1;
                        jN = 0;
                    end
                    nme_x_ = [  x_(i1:iN, :);
                                x_(nv+j1:nv+jN, :)  ];
                end

                %% call nonlinear function
                gg = cell(1, nargout);
                [gg{:}] = nme.(fcn)(nme_x_, sysx, nme_idx);

                %% insert the results in aggregate output args
                if sel
                    g(apidx, :) = gg{1};
                    if nargout > 1
                        if sysx
                            gv1(apidx, :) = gg{2};
                            gv2(apidx, :) = gg{3};
                            if nargout > 3 && nme.nz
                                gzr(apidx, :) = gg{4};
                                gzi(apidx, :) = gg{5};
                            end
                        else
                            gv1(apidx, i1:iN) = gg{2};
                            gv2(apidx, i1:iN) = gg{3};
                            if nargout > 3 && nme.nz
                                gzr(apidx, j1:jN) = gg{4};
                                gzi(apidx, j1:jN) = gg{5};
                            end
                        end
                    end
                else
                    g(i1:iN, :) = gg{1};
                    if nargout > 1
                        if sysx
                            gv1(i1:iN, :) = gg{2};
                            gv2(i1:iN, :) = gg{3};
                            if nargout > 3 && nme.nz
                                gzr(i1:iN, :) = gg{4};
                                gzi(i1:iN, :) = gg{5};
                            end
                        else
                            gv1(i1:iN, i1:iN) = gg{2};
                            gv2(i1:iN, i1:iN) = gg{3};
                            if nargout > 3 && nme.nz
                                gzr(i1:iN, j1:jN) = gg{4};
                                gzi(i1:iN, j1:jN) = gg{5};
                            end
                        end
                    end
                end
            end     %% for loop
        end

        function H = port_inj_nln_hess(obj, si, x_, lam, sysx, idx)
            % H = obj.port_inj_power_hess(x_, lam)
            % H = obj.port_inj_power_hess(x_, lam, sysx)
            % H = obj.port_inj_power_hess(x_, lam, sysx, idx)
            if nargin < 6
                idx = [];
                if nargin < 5
                    sysx = 1;
                end
            end

            %% current or power
            fcn = [si 'nln_hess'];
            fcn_list = [fcn '_list'];

            %% initialize
            n = 2 * length(x_);
            H = sparse(n, n);

            %% loop through elements w/gen nonlin Hessians, evaluate them
            pp = obj.get_idx('port');
            if ~sysx
                ss = obj.get_idx('state');
            end
            for kk = obj.(fcn_list)
                k = kk{1};      %% index into obj.elements
                nme = obj.elements{k};
                i1 = pp.i1.(nme.name)(1);
                iN = pp.iN.(nme.name)(end);

                %% set up x_ for nme & corresp row/col indices for nme
                if sysx
                    nme_x_ = x_;
                else
                    nv = obj.get_nv_(sysx);
                    nz = obj.nz;
                    if isfield(ss.i1, nme.name)
                        j1 = ss.i1.(nme.name)(1);
                        jN = ss.iN.(nme.name)(end);
                    else
                        j1 = 1;
                        jN = 0;
                    end
                    nme_x_ = [  x_(i1:iN, :);
                                x_(nv+j1:nv+jN, :)  ];

                    %% indices of rows/cols of H corresponding to nme x_
                    h = [(i1:iN) nv+(i1:iN) 2*nv+(j1:jN) 2*nv+nz+(j1:jN)].';
                end
                
                %% set up port index and lambda vectors for nme
                if ~isempty(idx)    %% selected ports only
                    apidx = find(idx >= i1 & idx <= iN);    %% aggregate port indices in range
                    if isempty(apidx)  %% skip if selected ports, but none in range
                        continue;
                    end
                    nme_idx = idx(apidx) - i1 + 1;  %% port index vector for nme
                    nme_lam = lam(apidx);           %% corresponding lam
                else                %% all ports
                    nme_idx = [];
                    nme_lam = lam(i1:iN);
                end

                %% call nonlinear function
                nme_H = nme.(fcn)(nme_x_, nme_lam, sysx, nme_idx);

                %% accumulate output
                if sysx
                    H = H + nme_H;
                else
                    H(h,h) = H(h,h) + nme_H;
                end
            end
        end

        function [G, Gv1, Gv2, Gzr, Gzi] = nodal_complex_current_balance(obj, x_)
            %% node incidence matrix
            C = obj.C;

            %% get port current injections with derivatives
            if nargout > 1
                [I, Iv1, Iv2, Izr, Izi] = obj.port_inj_current(x_, 1);
                Gv1 = C * Iv1;      %% Gva or Gvr
                Gv2 = C * Iv2;      %% Gvm or Gvi
                Gzr = C * Izr;
                Gzi = C * Izi;
            else
                I = obj.port_inj_current(x_, 1);
            end

            %% nodal current balance
            G = C * I;
        end

        function [G, Gv1, Gv2, Gzr, Gzi] = nodal_complex_power_balance(obj, x_)
            %% node incidence matrix
            C = obj.C;

            %% get port power injections with derivatives
            if nargout > 1
                [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1);
                Gv1 = C * Sv1;      %% Gva or Gvr
                Gv2 = C * Sv2;      %% Gvm or Gvi
                Gzr = C * Szr;
                Gzi = C * Szi;
            else
                S = obj.port_inj_power(x_, 1);
            end

            %% nodal power balance
            G = C * S;
        end

        function d2G = nodal_complex_current_balance_hess(obj, x_, lam)
            %% get port power injection hessians
            d2G = obj.port_inj_current_hess(x_, obj.C' * lam);
        end

        function d2G = nodal_complex_power_balance_hess(obj, x_, lam)
            %% get port power injection hessians
            d2G = obj.port_inj_power_hess(x_, obj.C' * lam);
        end

        function obj = port_inj_soln(obj)
            %% compute port injections
%             obj.soln.gi_ = obj.port_inj_current(obj.soln.x);
            obj.soln.gs_ = obj.port_inj_power(obj.soln.x);
        end

        function va = get_va(obj, idx)
            if isfield(obj.soln, 'v')           %% solved value
                if nargin < 2 || isempty(idx)
                    va = angle(obj.soln.v);
                else
                    va = angle(obj.soln.v(idx));
                end
            else                                %% initial value
                if nargin < 2 || isempty(idx)
                    va = obj.initial_voltage_angle();
                else
                    va = obj.initial_voltage_angle(idx);
                end
            end
        end

        %%-----  PF methods  -----
        function z_ = pf_update_z(obj, v_, z_, ad, Sinj, idx)
            %% update/allocate slack node active power injections
            %% and slack/PV node reactive power injections

            rpv = [ad.ref; ad.pv];      %% slack and PV nodes
            if nargin < 5 || isempty(Sinj)
                %% compute power injection at slack/PV nodes
                idx = find(any(obj.C(rpv, :), 1));  %% ports connected to slack/PV nodes
                Sinj = obj.port_inj_power([v_; z_], 1, idx);
            end
            Sref = obj.C(ad.ref, idx) * Sinj;
            Spv  = obj.C(ad.pv,  idx) * Sinj;
            Qpv = imag(Spv);

            %%-----  active power at slack nodes  -----
            %% coefficient matrix for power injection states
            CC = obj.C * obj.get_params([], 'N') * obj.D';

            %% coefficient matrix for power injection states for slack nodes
            CCref = CC(ad.ref, :);
            jr = find(any(CCref, 1));   %% indices of corresponding states

            %% active power injections at slack nodes
            Pref = real(Sref);

            %% allocate active power at slack nodes to 1st direct inj state
            %% find all z (except first one) with direct injection at each
            %% slack node
            [i, j] = find(CCref);
            if size(i, 2) > 1, i = i'; j = j'; end
            ij = sortrows([i j]);       %% 1st state comes 1st for each node
            [~, k1] = unique(ij(:, 1), 'first');%% index of 1st entry for each node
            %% all included states that are not 1st at their node
            jn = unique(ij(~ismember(1:length(i), k1), 2));

            %% if we have extra states (more than 1) for any node(s)
            if ~isempty(jn)
                %% augment update equation CC * (z - zprev) = -Pref with
                %% additional rows to force these states to remain fixed
                I = speye(obj.nz);
                CCref = [CCref; I(jn, :)];
                Pref = [Pref; zeros(length(jn), 1)];
            end

            %% update z for active injections at slack nodes
            z_(jr) = z_(jr) - CCref(:, jr) \ Pref;

            %%-----  reactive power at slack/PV nodes  -----
            %% coefficient matrix for power injection states for slack/PV nodes
            CCrpv = CC(rpv, :);
            jrpv = find(any(CCrpv, 1));
            Qrpv = [imag(Sref); Qpv] - CCrpv * imag(z_);

            %% find all z with direct injection at each slack/PV node
            [i, j] = find(CCrpv);
            if size(i, 2) > 1, i = i'; j = j'; end
            ij = sortrows([i j]);       %% 1st state comes 1st for each node
            [~, k1] = unique(ij(:, 1), 'first');%% index of 1st entry for each node
            % j1 = ij(k1, 2);     %% indices of states that are 1st at their node
            kn = find(~ismember(1:length(i), k1));  %% indices of entries that are not first state for corresponding node
            %% all included states that are not 1st at their node
            jn = unique(ij(kn, 2));
            in = unique(ij(kn, 1));     %% nodes with multiple states

            %% if we have extra states (more than 1) for some node(s)
            if ~isempty(jn)
                %% find ranges for relevant state vars to allocate reactive
                %% power proportional to ranges
                [~, mn, mx] = obj.params_var('zi');

                %% This code is currently not designed to handle injections
                %% from states that affect more than a single node, so we
                %% throw a warning if we find such a case
                if any(sum(CCrpv ~= 0) > 1)
                    k = find(sum(CCrpv ~= 0) > 1);
                    warning('pf_update_z:multiple_nodes', ...
                        'mp_network_ac/pf_update_z: unable to distribute reactive power due to z var %d affecting multiple nodes.', k(1));
                end

                %% define a numerical proxy to replace +/- Inf limits
                %% start by setting M equal to average injection at node
                %% CCrpv' * Qrpv is total Qrpv at each corresponding state
                %% CCrpv' * sum(CCrpv, 2) is the number of injections at node
                M = abs((CCrpv' * Qrpv) ./ (CCrpv' * sum(CCrpv, 2)));
                %% add abs value of upper and lower bound to avg nodal injection
                M(~isinf(mx)) = M(~isinf(mx)) + abs(mx(~isinf(mx)));
                M(~isinf(mn)) = M(~isinf(mn)) + abs(mn(~isinf(mn)));
                %% set M for each state to sum over all states at same node
                M = CC' * CC * M;
                %% replace +/- Inf limits with proxy +/- M
                mn(mn ==  Inf) =  M(mn ==  Inf);
                mn(mn == -Inf) = -M(mn == -Inf);
                mx(mx ==  Inf) =  M(mx ==  Inf);
                mx(mx == -Inf) = -M(mx == -Inf);

                %% find indices of states with largest range at its node
                r = mx - mn;    %% size of range for each state
                [rmax, j] = max(-CCrpv * spdiags(r, 0, obj.nz, obj.nz), [], 2);

                %% set ranges to 1 for states at nodes where all ranges are 0
                %% (results in equal limit violations)
                %% find nodes where all corresponding ranges are zero
                i0 = find(abs(rmax) < 10*eps);
                if ~isempty(i0)
                    rmax(i0) = 1;           %% set these ranges to one
                    j0 = find(any(CCrpv(i0, :), 1));    %% corresponding states
                    r(j0) = -CCrpv(i0, j0)' * rmax(i0); %% apply to all states at same node
                end

                %% augment update equation ...
                %%  CCrpv * (z - zprev + 1j * imag(zprev)) = -1j * Qrpv
                %% with additional rows to force these states to allocate in
                %% proportion to min-max range, i.e. all states k at node have
                %%  q_k = -qmin_k - r_k * lam, where lam is same for all k at node
                %% we solve for lam using eqn with largest r_k, then
                %% substitute in other equations to get the set of constraints
                %% to add
                R = sparse(obj.nz, obj.nz);
                b = zeros(obj.nz, 1);
                jn = [];    %% initialize list of states that do not have max range at their node
                for i = 1:length(in)    %% for all nodes with multiple states
                    ii = in(i);         %% node ii, with multiple states
                    jj = ij(ij(:, 1) == ii, 2); %% states at node ii
                    rr = r(jj) / r(j(ii));  %% range of each z / max z range
                    R(jj, j(ii)) = rr;  %% rows for states @ node ii, col of state w/max range
                    b(jj) = b(jj) + rr * mn(j(ii)) - mn(jj);
                    jj(jj == j(ii)) = [];
                    jn = [jn; jj];  %% add states at node i that do not have max range
                end
                A = speye(obj.nz) - R;
                Qrpv = [Qrpv; b(jn)];
                CCrpv = [CCrpv; A(jn, :)];
            end

            %% update z for reactive injections at slack/PV nodes
            z0 = z_(jrpv);
            z_(jrpv) = z0 - 1j * (CCrpv(:, jrpv) \ Qrpv + imag(z0));
        end

        %%-----  CPF methods  -----
        function [vx_, z_, x_] = cpf_convert_x(obj, mmx, ad, only_v)
            nmt = ad.nmt;
            lam = mmx(end);     %% continuation parameter lambda

            %% update voltages and get base z_
            [vx_,  zb_] = obj.pf_convert_x(mmx(1:end-1), ad,     1);
            [vxt_, zt_] = nmt.pf_convert_x(mmx(1:end-1), ad.adt, 1);
            assert(norm(vx_-vxt_, Inf) < eps);

            %% compute z_ as function of continuation parameter lambda
            z_ = (1-lam) * zb_+ lam * zt_;

            %% update dependent portions of z, if requested
            if nargin < 4 || ~only_v
                rpv = [ad.ref; ad.pv];      %% slack and PV nodes
                idx = find(any(obj.C(rpv, :), 1));  %% ports connected to slack/PV nodes
                Sinjb = obj.port_inj_power([vx_; zb_], 1, idx);
                Sinjt = nmt.port_inj_power([vx_; zt_], 1, idx);
                Sinj = (1-lam) * Sinjb + lam * Sinjt;
                z_ = obj.pf_update_z(vx_, z_, ad, Sinj, idx);
            end

            if nargout < 2
                vx_ = [vx_; z_];
            elseif nargout > 2
                x_ = [vx_; z_];
            end
        end

        function ad = cpf_check_xfer(obj, nmt, ad, adt)
            %% ensure base and target parameters are identical,
            %% except for fixed power injections & states
            [L,  N,  i,  s ] = obj.get_params([], {'L', 'N', 'i', 's'});
            [Lt, Nt, it, st] = nmt.get_params([], {'L', 'N', 'i', 's'});
            tol = 1e-12;

            %% nodal transfer from direct power injection states
            NN = obj.C * N;
            NNt = obj.C * Nt;
            LL = obj.C * L;
            LLt = obj.C * Lt;
            zzr = ad.zr - adt.zr;
            zzi = ad.zi - adt.zi;
            kref = find(any(LL(ad.ref, :), 1) | any(NN(ad.ref, :), 1));
            kpv  = find(any(LL(ad.pv,  :), 1) | any(NN(ad.pv,  :), 1));
            zzr(kref) = 0;   %% zero active transfer at slack node
            zzi(kpv)  = 0;   %% zero reactive transfer at PV nodes
            zz = zzr + 1j * zzi;

            %% nodal transfer from constant power injections
            ss = obj.C * (s - st);

            %% create transfer vector from diff in direct power injections
            %% between base and target, from constant power elements and
            %% direct power injection states, used only to auto select
            %% largest transfer bus for default voltage nose-curve plot
            ad.xfer = ss + NN * zz;

            %% Power flow equations must be linear in continuation parameter
            %% lambda. To do that we must ensure that ...
            %% 1. Specified voltages do not vary with lambda.
            [va,  vm ] = obj.aux_data_va_vm(ad);
            [vat, vmt] = obj.aux_data_va_vm(adt);
            rpv = [ad.ref; ad.pv];
            if norm(va(ad.ref)-vat(ad.ref), Inf) > tol
                error('mpe_network_ac/cpf_check_xfer: base and target cases must have identical voltages angles at reference nodes.')
            end
            if norm(vm(rpv)-vmt(rpv), Inf) > tol
                error('mpe_network_ac/cpf_check_xfer: base and target cases must have identical voltage magnitudes at reference and PV nodes.')
            end
            %% 2. Elements of z that vary with lambda must have only constant
            %%    coefficients, i.e. corresponding columns of L and N must be
            %%    identical in base and target models.
            k = find(zz);
            if norm(LL(:,k) - LLt(:,k), Inf) > tol
                error('mpe_network_ac/cpf_check_xfer: base and target cases must have identical coefficients for any current injection state variables that vary from base to target.')
            end
            if norm(NN(:,k) - NNt(:,k), Inf) > tol
                error('mpe_network_ac/cpf_check_xfer: base and target cases must have identical coefficients for any power injection state variables that vary from base to target.')
            end
        end

        function obj = cpf_data_model_update(obj, mm, nm, dm, mpopt)
            %% each element updates its data model
            for k = 1:length(obj.elements)
                obj.elements{k}.cpf_data_model_update(mm, nm, dm, mpopt);
            end
        end

        function [names, vals] = cpf_pne_output_fcn(obj, ad, x, x_hat)
            %% [names, vals] = obj.cpf_pne_history(ad, x, x_hat)
            %% names = obj.cpf_pne_history(ad)
            names = {'V_hat', 'V'};
            if nargin > 2
                [V_hat, ~] = obj.cpf_convert_x(x_hat, ad, 1);
                [V,     ~] = obj.cpf_convert_x(x,     ad, 1);
                vals = {V_hat, V};
            end
        end

        function y = cpf_plot_yfcn(obj, dm, ad, v_, bus_num)
            %% find node idx from external bus number
            b2i = dm.elements.bus.ID2i;     %% bus num to idx mapping
            nidx = obj.get_node_idx('bus');

            k = find( bus_num < 0 | bus_num > length(b2i) );
            if ~isempty(k)
                error('mpe_network/cpf_plot_yfcn: %d is not a valid bus number for CPF voltage plot', bus_idx(k));
            end

            idx = nidx(b2i(bus_num));
            y = abs(v_(idx, :));
        end

        function idx = cpf_plot_idx_default(obj, dm, ad)
            %% plot voltage of PQ node with max transfer as default
            nidx = obj.get_node_idx('bus');     %% node indices of buses
            [~, i] = max(abs(ad.xfer(ad.pq)) .* ismember(ad.pq, nidx));
            bi = ad.pq(i);                      %% index of bus w/max transfer
            idx = dm.elements.bus.ID(bi);       %% bus num of same bus
        end

        function opt = cpf_add_callbacks(obj, opt, mm, dm, mpopt)
            qlim = mpopt.cpf.enforce_q_lims;    %% enforce reactive limits
            plim = mpopt.cpf.enforce_p_lims;    %% enforce active limits
            vlim = mpopt.cpf.enforce_v_lims;    %% enforce voltage magnitude limits
            flim = mpopt.cpf.enforce_flow_lims; %% enforce branch flow limits
            
            %% initialize event and callback options
            if ~isfield(opt, 'events') || isempty(opt.events)
                opt.events = {};
            end
            if ~isfield(opt, 'callbacks') || isempty(opt.callbacks)
                opt.callbacks = {};
            end

            if flim
                opt.events{end+1} = { ...
                    'FLIM', ...
                    @(cx, opt)cpf_event_flim(obj, cx, opt, mm, dm, mpopt), ...
                    mpopt.cpf.flow_lims_tol };
                opt.callbacks{end+1} = { ...
                    @(k, nx, cx, px, s, opt)cpf_callback_flim(obj, k, nx, cx, px, s, opt, mm, dm, mpopt), ...
                    53 };
            end

            if vlim
                opt.events{end+1} = { ...
                    'VLIM', ...
                    @(cx, opt)cpf_event_vlim(obj, cx, opt, mm, dm, mpopt), ...
                    mpopt.cpf.v_lims_tol };
                opt.callbacks{end+1} = { ...
                    @(k, nx, cx, px, s, opt)cpf_callback_vlim(obj, k, nx, cx, px, s, opt, mm, dm, mpopt), ...
                    52 };
            end

            if qlim
                opt.events{end+1} = { ...
                    'QLIM', ...
                    @(cx, opt)cpf_event_qlim(obj, cx, opt, mm, dm, mpopt), ...
                    mpopt.cpf.q_lims_tol };
                opt.callbacks{end+1} = { ...
                    @(k, nx, cx, px, s, opt)cpf_callback_qlim(obj, k, nx, cx, px, s, opt, mm, dm, mpopt), ...
                    41 };
            end

            if plim
                opt.events{end+1} = { ...
                    'PLIM', ...
                    @(cx, opt)cpf_event_plim(obj, cx, opt, mm, dm, mpopt), ...
                    mpopt.cpf.p_lims_tol };
                opt.callbacks{end+1} = { ...
                    @(k, nx, cx, px, s, opt)cpf_callback_plim(obj, k, nx, cx, px, s, opt, mm, dm, mpopt), ...
                    40 };
            end
        end

        function efv = cpf_event_flim(obj, cx, opt, mm, dm, mpopt)
            %% get branch flow constraints
            branch_nme = obj.elements.branch;
            branch_dme = branch_nme.data_model_element(dm);
            rate_a = branch_dme.rate_a * dm.base_mva;
            ibr = find(rate_a ~= 0 & rate_a < 1e10);
            nl2 = length(ibr);          %% number of constrained branches

            if nl2
                nl = branch_nme.nk;     %% port indexes

                %% convert cx.x back to x_
                ad = mm.aux_data;
                x_ = obj.cpf_convert_x(cx.x, ad);

                %% branch flows
                S_fr = branch_nme.port_inj_power(x_, 1, ibr)    * dm.base_mva;
                S_to = branch_nme.port_inj_power(x_, 1, nl+ibr) * dm.base_mva;
                S_fr = sqrt(S_fr .* conj(S_fr));
                S_to = sqrt(S_to .* conj(S_to));

                %% branch flow lim event function
                efv = max(S_fr, S_to) - rate_a(ibr);
            else
                efv = NaN;
            end
        end

        function efv = cpf_event_qlim(obj, cx, opt, mm, dm, mpopt)
            ad = mm.aux_data;

            %% convert cx.x back to v_, z_
            [v_, z_] = obj.cpf_convert_x(cx.x, ad);

            %% coefficient matrix for power injection states
            rpv = [ad.ref; ad.pv];      %% slack and PV nodes
            CCrpv = obj.C(rpv, :) * obj.get_params([], 'N') * obj.D';
            jrpv = find(any(CCrpv, 1)); %% indices for states at slack/PV nodes

            %% limit violations at slack/PV nodes
            [~, zi_min, zi_max] = obj.params_var('zi'); %% bounds on zi
            v_Qmax = NaN(obj.nz, 1);
            v_Qmin = v_Qmax;
            v_Qmax(jrpv) = imag(z_(jrpv)) - zi_max(jrpv);
            v_Qmin(jrpv) = zi_min(jrpv) - imag(z_(jrpv));

            %% assemble event function value
            efv = [v_Qmax; v_Qmin] * dm.base_mva;
        end

        function efv = cpf_event_plim(obj, cx, opt, mm, dm, mpopt)
            ad = mm.aux_data;

            %% convert cx.x back to v_, z_
            [v_, z_] = obj.cpf_convert_x(cx.x, ad);

            %% limit violations
            [~, ~, zr_max] = obj.params_var('zr'); %% bounds on zr
            v_Pmax = real(z_) - zr_max;

            %% ignore those that are already at their max limit
            if isfield(cx.cbs, 'plim') && ~isempty(cx.cbs.plim.idx)
                v_Pmax(cx.cbs.plim.idx) = NaN;
            end

            %% assemble event function value
            efv = v_Pmax * dm.base_mva;
        end

        function [nx, cx, s] = cpf_callback_flim(obj, k, nx, cx, px, s, opt, mm, dm, mpopt)
            %% initialize
            if k == 0   %% check for base case flow violations
                %% get branch flow constraints
                branch_nme = obj.elements.branch;
                branch_dme = branch_nme.data_model_element(dm);
                rate_a = branch_dme.rate_a * dm.base_mva;
                ibr = find(rate_a ~= 0 & rate_a < 1e10);
                nl2 = length(ibr);          %% number of constrained branches

                if nl2
                    nl = branch_nme.nk;     %% port indexes

                    %% convert cx.x back to x_
                    ad = mm.aux_data;
                    x_ = obj.cpf_convert_x(cx.x, ad);

                    %% branch flows
                    S_fr = branch_nme.port_inj_power(x_, 1, ibr)    * dm.base_mva;
                    S_to = branch_nme.port_inj_power(x_, 1, nl+ibr) * dm.base_mva;
                    S_fr = sqrt(S_fr .* conj(S_fr));
                    S_to = sqrt(S_to .* conj(S_to));

                    %% violated branch flows
                    if any(max(S_fr, S_to) > rate_a(ibr))
                        %% find the lines and which lim(s)
                        iL = find(max(S_fr, S_to) > rate_a(ibr));
                        msg = '';
                        for j = 1:length(iL)
                            L = ibr(iL(j));
                            fidx = find(branch_nme.C(:, L));
                            tidx = find(branch_nme.C(:, nl+L));
                            flabel = obj.set_type_label('node', fidx, dm);
                            tlabel = obj.set_type_label('node', tidx, dm);

                            msg = sprintf('%sbranch flow limit violated in base case: %s -- %s exceeds limit of %g MVA\n',...
                                msg, flabel, tlabel, rate_a(L));
                        end

                        %% prepare to terminate
                        s.done = 1;
                        s.done_msg = msg;
                    end
                end
            end

            %% skip if finalize or done
            if k < 0 || s.done
                return;
            end

            %% handle event
            ev = pne_detected_event(s.events, 'FLIM', 1);   %% zero only
            if ~isempty(ev)
                if opt.verbose > 3
                    msg = sprintf('%s\n    ', ev.msg);
                else
                    msg = '';
                end

                %% get branch flow constraints
                branch_nme = obj.elements.branch;
                branch_dme = branch_nme.data_model_element(dm);
                rate_a = branch_dme.rate_a * dm.base_mva;
                ibr = find(rate_a ~= 0 & rate_a < 1e10);
                nl2 = length(ibr);      %% number of constrained branches
                nl = branch_nme.nk;     %% port indexes

                %% find branch(es) with violated lim(s)
                iL = ev.idx;            %% event function index
                for j = 1:length(iL)
                    L = ibr(iL(j)); %% index of critical branch event of interest
                    fidx = find(branch_nme.C(:, L));
                    tidx = find(branch_nme.C(:, nl+L));
                    flabel = obj.set_type_label('node', fidx, dm);
                    tlabel = obj.set_type_label('node', tidx, dm);

                    msg = sprintf('%sbranch flow limit reached\nbranch: %s -- %s at limit of %g MVA @ lambda = %.4g, in %d continuation steps',...
                        msg, flabel, tlabel, rate_a(L), nx.x(end), k);
                end

                %% prepare to terminate
                s.done = 1;
                s.done_msg = msg;
            end
        end

        function [nx, cx, s] = cpf_callback_qlim(obj, k, nx, cx, px, s, opt, mm, dm, mpopt)
            %% skip if initialize, finalize or done
            if k <= 0 || s.done
                return;
            end

            %% handle event
            [ev, i] = pne_detected_event(s.events, 'QLIM', 1);  %% zero only
            if ~isempty(ev)
                ad = mm.aux_data;
                if ad.nref ~= 1
                    error('mp_network_acps/cpf_callback_qlim: ''cpf.enforce_qlims'' option only valid for systems with exactly one REF bus');
                end

                efidx = ev.idx;             %% event function index
                [~, zi_min, zi_max] = obj.params_var('zi'); %% bounds on zi
                if opt.verbose > 3
                    msg = sprintf('%s\n    ', ev.msg);
                else
                    msg = '';
                end
                for j = 1:length(efidx)
                    %% find index of z var
                    idx = efidx(j);         %% index of z var
                    if idx <= obj.nz
                        maxlim = 1;         %% Qmax violation
                        lim = zi_max(idx) * dm.base_mva;
                        lim_type = 'Qmax';
                    else
                        idx = idx - obj.nz; %% correct index of z var
                        maxlim = 0;         %% Qmin violation
                        lim = zi_min(idx) * dm.base_mva;
                        lim_type = 'Qmin';
                    end

                    %% get label for z var
                    zlabel = obj.set_type_label('state', idx, dm);

                    %% get label for corresponding node
                    CC = obj.C * obj.get_params([], 'N') * obj.D';
                    nidx = find(CC(:, idx));
                    nlabel = obj.set_type_label('node', nidx, dm);

                    msg = sprintf('%s%s @ %s reached %g MVAr %s lim @ lambda = %.4g : %s converted to PQ', ...
                            msg, zlabel, nlabel, lim, lim_type, ...
                            nx.x(end), nlabel);

                    %% set Q to exact limit
                    [v_, z_] = obj.cpf_convert_x(nx.x, ad);
                    z_(idx) = real(z_(idx)) + 1j * lim / dm.base_mva;

                    %% change node type to PQ
                    obj.set_node_type_pq(dm, nidx);

                    %% check for existence of remaining slack/PV bus
                    try
                        %% potentially pick new reference bus
                        [ref, pv, pq] = obj.node_types(obj, dm);
                    catch
                        s.done = 1;
                        s.done_msg = 'No REF or PV nodes remaining.';

                        %% undo change of last REF to PQ
                        obj.set_node_type_ref(dm, ad.ref);

                        break;
                    end

                    if ~s.done
                        %% get target case data, network models
                        dmt = ad.dmt;
                        nmt = ad.nmt;

                        %% change node type in target case
                        nmt.set_node_type_pq(dmt, nidx);

                        %% zero out Q transfer for bus
                        ss = obj.set_type_idx_map('zi', idx);
                        obj.zi.data.v0.(ss.name)(ss.i) = imag(z_(idx));
                        nmt.zi.data.v0.(ss.name)(ss.i) = imag(z_(idx));

                        %% if slack changed ...
                        if ref ~= ad.ref
                            %% find zr corresponding to all zr at ref node
                            zref = find(CC(ad.ref, :));
                            ss = obj.set_type_idx_map('zr', zref);

                            %% zero out P transfer at old ref
                            for kk = 1:length(ss)
                                obj.zr.data.v0.(ss(kk).name)(ss(kk).i) = real(z_(zref(kk)));
                                nmt.zr.data.v0.(ss(kk).name)(ss(kk).i) = real(z_(zref(kk)));
                            end

                            %% update voltage angle at new ref node
                            ss = obj.set_type_idx_map('va', ref);
                            obj.va.data.v0.(ss.name)(ss.i) = angle(v_(ref));
                            nmt.va.data.v0.(ss.name)(ss.i) = angle(v_(ref));
                        end
                    end
                end
                if ~s.done
                    s.done = 1;
                    dir_from_jac_eigs = isempty(find(nx.z(ad.npv+ad.npq+1:end-1) > 0, 1));
                    s.warmstart = struct('nmt', nmt, 'dmt', dmt, ...
                        'dir_from_jac_eigs', dir_from_jac_eigs);
                end
                s.events(i).msg = msg;
            end
        end

        function [nx, cx, s] = cpf_callback_plim(obj, k, nx, cx, px, s, opt, mm, dm, mpopt)
            %% skip if finalize or done
            if k < 0 || s.done
                return;
            elseif k == 0 && ~isfield(cx.cbs, 'plim')
                cx.cbs.plim.idx = [];
            end

            %% handle event
            [ev, i] = pne_detected_event(s.events, 'PLIM', 1);  %% zero only
            if ~isempty(ev)
                ad = mm.aux_data;
                if ad.nref ~= 1
                    error('mp_network_acps/cpf_callback_plim: ''cpf.enforce_plims'' option only valid for systems with exactly one REF bus');
                end

                efidx = ev.idx;             %% event function index
                [~, ~, zr_max] = obj.params_var('zr'); %% bounds on zr
                if opt.verbose > 3
                    msg = sprintf('%s\n    ', ev.msg);
                else
                    msg = '';
                end
                for j = 1:length(efidx)
                    %% find index of z var
                    idx = efidx(j);         %% index of z var
                    lim = zr_max(idx) * dm.base_mva;

                    %% get label for z var
                    zlabel = obj.set_type_label('state', idx, dm);

                    %% get label for corresponding node
                    CC = obj.C * obj.get_params([], 'N') * obj.D';
                    nidx = find(CC(:, idx));
                    nlabel = obj.set_type_label('node', nidx, dm);

                    msg = sprintf('%s%s @ %s reached %g MW pg upper bound @ lambda = %.4g', ...
                            msg, zlabel, nlabel, lim, nx.x(end));

                    %% set P to exact limit
                    [v_, z_] = obj.cpf_convert_x(nx.x, ad);
                    z_(idx) = lim / dm.base_mva + 1j * imag(z_(idx));

                    dmt = ad.dmt;
                    nmt = ad.nmt;

                    %% find zr corresponding to all zr at this node
                    k = find(CC(nidx, :));
                    ss = obj.set_type_idx_map('zr', k);

                    %% set z to limit at this node
                    for kk = 1:length(ss)
                        obj.zr.data.v0.(ss(kk).name)(ss(kk).i) = real(z_(k(kk)));
                        nmt.zr.data.v0.(ss(kk).name)(ss(kk).i) = real(z_(k(kk)));
                    end

                    %% save index of element at limit
                    nx.cbs.plim.idx = [nx.cbs.plim.idx; idx];

                    %% if it is at the ref node
                    if nidx == ad.ref
                        %% find free injections (not at max lim)
                        free = ones(obj.nz, 1);
                        free(nx.cbs.plim.idx) = 0;

                        %% first node with any free injection
                        ref = find(any(CC * spdiags(free, 0, obj.nz, obj.nz), 2));
                        if isempty(ref)
                            s.done = 1;
                            s.done_msg = 'All generators at Pmax';
                            break;
                        else
                            %% convert this node to PV, new ref bus to REF
                            obj.set_node_type_pv(dm, nidx);
                            nmt.set_node_type_pv(dmt, nidx);
                            obj.set_node_type_ref(dm, ref);
                            nmt.set_node_type_ref(dmt, ref);

                            %% update voltage angle at new ref node
                            ss = obj.set_type_idx_map('va', ref);
                            obj.va.data.v0.(ss.name)(ss.i) = angle(v_(ref));
                            nmt.va.data.v0.(ss.name)(ss.i) = angle(v_(ref));

                            rlabel = obj.set_type_label('node', ref, dm);
                            msg = sprintf('%s : ref changed from %s to %s', ...
                                msg, nlabel, rlabel);
                        end
                    end
                end
                if ~s.done
                    s.done = 1;
                    s.warmstart = struct('nmt', nmt, 'dmt', dmt);
                end
                s.events(i).msg = msg;
            end
        end
    end     %% methods
end         %% classdef
