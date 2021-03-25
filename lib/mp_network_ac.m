classdef mp_network_ac < mp_network% & mp_form_ac
%MP_NETWORK_AC Abstract class, explicitly a subclass of MP_NETWORK and
%              implicitly assumed to be subclasses of MP_FORM_AC as well

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
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
            for k = 1:length(obj.elm_list)
                if ~isempty(obj.elm_list{k}.inln)
                    obj.inln_list{end+1} = k;
                    if ~isempty(obj.elm_list{k}.inln_hess)
                        obj.inln_hess_list{end+1} = k;
                    end
                end
                if ~isempty(obj.elm_list{k}.snln)
                    obj.snln_list{end+1} = k;
                    if ~isempty(obj.elm_list{k}.snln_hess)
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
                k = kk{1};      %% index into obj.elm_list
                nme = obj.elm_list{k};
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
                k = kk{1};      %% index into obj.elm_list
                nme = obj.elm_list{k};
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


        %%-----  PF methods  -----
        function ad = pf_aux_data(obj, dm, mpopt)
            %% get model variables
            vvars = obj.model_vvars();
            zvars = obj.model_zvars();
            v1 = obj.params_var(vvars{1});
            v2 = obj.params_var(vvars{2});
            zr = obj.params_var(zvars{1});
            zi = obj.params_var(zvars{2});

            %% get node types
            [ref, pv, pq] = obj.node_types(obj, dm);

            %% create aux_data struct
            ad = struct( ...
                'v1', v1, ...               %% initial value of v1 (va or vr)
                'v2', v2, ...               %% initial value of v2 (vm or vi)
                'zr', zr, ...               %% initial value of zr
                'zi', zi, ...               %% initial value of zi
                'ref',  ref, ...            %% REF node indices
                'nref', length(ref), ...    %% number of REF nodes
                'pv',  pv, ...              %% PV node indices
                'npv', length(pv), ...      %% number of PV nodes
                'pq',  pq, ...              %% PQ node indices
                'npq', length(pq) ...       %% number of PQ nodes
            );
        end

        function opt = pf_solve_opts(obj, mm, dm, mpopt)
            switch mpopt.pf.alg
                case 'DEFAULT'
                    opt = mpopt2nleqopt(mpopt, mm.problem_type(), 'DEFAULT');
                case {'NR', 'NR-SP', 'NR-SC', 'NR-SH', 'NR-IP', 'NR-IC', 'NR-IH'}
                    opt = mpopt2nleqopt(mpopt, mm.problem_type(), 'NEWTON');
                case {'FDXB', 'FDBX'}
                    opt = mpopt2nleqopt(mpopt, mm.problem_type(), 'FD');
                    opt.fd_opt.jac_approx_fcn = @()obj.pf_fd_jac_approx(mm, dm, mpopt);
                    opt.fd_opt.labels = {'P', 'Q'};
                case 'FSOLVE'
                    opt = mpopt2nleqopt(mpopt, mm.problem_type(), 'FSOLVE');
                case 'GS'
                    opt = mpopt2nleqopt(mpopt, mm.problem_type(), 'GS');
                    opt.gs_opt.x_update_fcn = ...
                        @(x, f)obj.pf_gs_x_update(x, f, mm, dm, mpopt);
                case 'ZG'
                    opt = mpopt2nleqopt(mpopt, mm.problem_type(), 'ZG');
                    zg_x_update = @(x, f)obj.pf_zg_x_update(x, f, mm, dm, mpopt);
                    opt.core_sp = struct(...
                        'alg',              'ZG', ...
                        'name',             'Implicit Z-bus Gauss', ...
                        'default_max_it',   1000, ...
                        'need_jac',         0, ...
                        'update_fcn',       zg_x_update  );
                otherwise
                    error('mp_network_ac/pf_solve_opts: invalid value for MPOPT.PF.ALG (%s)', mpopt.pf.alg);
            end
            opt.verbose = mpopt.verbose;
        end

        function z_ = pf_update_z(obj, x_, z_, ad)
            %% update/allocate slack node active power injections
            %% and slack/PV node reactive power injections

            %% coefficient matrix for power injection states
            CC = obj.C * obj.get_params([], 'N') * obj.D';

            %% compute power injection at slack/PV nodes
            rpv = [ad.ref; ad.pv];      %% slack and PV nodes
            idx = find(any(obj.C(rpv, :), 1));  %% ports connected to slack/PV nodes
            Sinj = obj.port_inj_power([x_; z_], 1, idx);
            Sref = obj.C(ad.ref, idx) * Sinj;
            Spv  = obj.C(ad.pv,  idx) * Sinj;

            %%-----  active power at slack nodes  -----
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
            Qrpv = imag([Sref; Spv]) - CCrpv * imag(z_);

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
                %% throw an error if we find such a case
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
        function ad = cpf_aux_data(obj, nmt, dm, dmt, mpopt)
            ad = obj.pf_aux_data(dm, mpopt);

            ad.xfer = obj.cpf_xfer(nmt);
            ad.nmt = nmt;
            ad.dmt = dmt;
            ad.mpopt = mpopt;
        end

        function xfer = cpf_xfer(obj, nmt)
            %% ensure base and target parameters are identical,
            %% except for fixed power injections & states
            [Y,  L,  M,  N,  i,  s ] = obj.get_params();
            [Yt, Lt, Mt, Nt, it, st] = nmt.get_params();
            tol = 1e-12;
            if norm(Y-Yt, Inf) > tol || norm(L-Lt, Inf) > tol || ...
                    norm(M-Mt, Inf) > tol || norm(N-Nt, Inf) > tol || ...
                    norm(i-it, Inf) > tol
                error('mpe_network_ac/cpf_aux_data: base and target cases must differ only in direct power injections')
            end

            %% create transfer vector from diff between base & target cases
            z0  = obj.params_var('zr') + 1j * obj.params_var('zi');
            z0t = nmt.params_var('zr') + 1j * nmt.params_var('zi');
            sg = N * (z0 - z0t);

            xfer = obj.C * (s - st + sg);
        end

        function cpf_add_vars(obj, mm, nm, dm, mpopt)
            obj.pf_add_vars(mm, nm, dm, mpopt);
            mm.add_var('lambda', 1, 0);
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
            
            if qlim
                opt.events{end+1} = { ...
                    'QLIM', ...
                    @(cx, opt)cpf_event_qlim(obj, cx, opt, mm, dm, mpopt), ...
                    mpopt.cpf.q_lims_tol };
                opt.callbacks{end+1} = { ...
                    @(k, nx, cx, px, s, opt)cpf_callback_qlim(obj, k, nx, cx, px, s, opt, mm, dm, mpopt), ...
                    41 };
            end
        end


        %%-----  OPF methods  -----
        function [g, dg] = opf_current_balance_fcn(obj, x_)
            if nargout > 1
                [G, Gv1, Gv2, Gzr, Gzi] = obj.nodal_complex_current_balance(x_);
                Gx = [Gv1 Gv2 Gzr Gzi];
                dg = [  real(Gx);       %% Re{I} mismatch w.r.t v1, v2, zr, zi
                        imag(Gx)    ];  %% Im{I} mismatch w.r.t v1, v2, zr, zi
            else
                G = obj.nodal_complex_current_balance(x_);
            end
            g = [ real(G);              %% real current mismatch
                  imag(G) ];            %% imaginary current mismatch
        end

        function [g, dg] = opf_power_balance_fcn(obj, x_)
            if nargout > 1
                [G, Gv1, Gv2, Gzr, Gzi] = obj.nodal_complex_power_balance(x_);
                Gx = [Gv1 Gv2 Gzr Gzi];
                dg = [  real(Gx);       %% P mismatch w.r.t v1, v2, zr, zi
                        imag(Gx)    ];  %% Q mismatch w.r.t v1, v2, zr, zi
            else
                G = obj.nodal_complex_power_balance(x_);
            end
            g = [ real(G);              %% active power (P) mismatch
                  imag(G) ];            %% reactive power (Q) mismatch
        end

        function d2G = opf_current_balance_hess(obj, x_, lam)
            nlam = length(lam) / 2;
            lamIr = lam(1:nlam);
            lamIi = lam((1:nlam)+nlam);

            d2Gr = obj.nodal_complex_current_balance_hess(x_, lamIr);
            d2Gi = obj.nodal_complex_current_balance_hess(x_, lamIi);

            d2G = real(d2Gr) + imag(d2Gi);
        end

        function d2G = opf_power_balance_hess(obj, x_, lam)
            nlam = length(lam) / 2;
            lamP = lam(1:nlam);
            lamQ = lam((1:nlam)+nlam);

            d2Gr = obj.nodal_complex_power_balance_hess(x_, lamP);
            d2Gi = obj.nodal_complex_power_balance_hess(x_, lamQ);

            d2G = real(d2Gr) + imag(d2Gi);
        end

        function opf_add_system_costs(obj, mm, dm, mpopt)
            %% can be overridden to add additional system costs

            %% legacy user-defined costs
            obj.opf_add_legacy_user_costs(mm, dm, 0);
        end

        function opf_add_legacy_user_constraints(obj, mm, dm, mpopt)
            %% call parent
            opf_add_legacy_user_constraints@mp_network(obj, mm, dm, mpopt);

            uc = dm.opf_legacy_user_constraints();
            for k = 1:length(uc)
                mm.add_nln_constraint(uc{k}{:});
            end
        end
    end     %% methods
end         %% classdef
