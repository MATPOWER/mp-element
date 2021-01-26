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

        function obj = pf_add_vars(obj, mm, nm, dm, mpopt)
            %% get model variables
            vvars = obj.model_vvars();

            %% index vectors
            ad = mm.get_userdata('aux_data');
            pvq = [ad.pv; ad.pq];

            %% voltage angles
            st = obj.(vvars{1});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    mm.add_var(name, ad.npv+ad.npq, d.v0.(name)(pvq), d.vl.(name)(pvq), d.vu.(name)(pvq));
                else
                    error('mp_network_acps/pf_add_vars: handling of indexed sets not implmented here (yet)');
                end
            end

            %% voltage magnitudes
            st = obj.(vvars{2});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    mm.add_var(name, ad.npq, d.v0.(name)(ad.pq), d.vl.(name)(ad.pq), d.vu.(name)(ad.pq));
                else
                    error('mp_network_acps/pf_add_vars: handling of indexed sets not implmented here (yet)');
                end
            end
        end

        function [vx_, z_, x_] = pf_convert_x(obj, mmx, ad, only_v)
            %% x = obj.pf_convert(mmx, ad)
            %% [v, z] = obj.pf_convert(mmx, ad)
            %% [v, z, x] = obj.pf_convert(mmx, ad)
            %% ... = obj.pf_convert(mmx, ad, only_v)

            %% update v_, z_ from mmx
            ad.v1([ad.pv; ad.pq]) = mmx(1:ad.npv+ad.npq);                   %% va
            ad.v2(ad.pq)          = mmx(ad.npv+ad.npq+1:ad.npv+2*ad.npq);   %% vm
            vx_ = ad.v2 .* exp(1j * ad.v1);
            z_ = ad.zr + 1j * ad.zi;

            %% update z, if requested
            if nargin < 4 || ~only_v
                z_ = obj.pf_update_z(vx_, z_, ad);
            end

            if nargout < 2
                vx_ = [vx_; z_];
            elseif nargout > 2
                x_ = [vx_; z_];
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
                [S, Sva, Svm] = obj.port_inj_power([v_; z_], 1);

                SSva = C * Sva;
                SSvm = C * Svm;
                J = [   real(SSva(pvq,   pvq))  real(SSvm(pvq,   ad.pq));
                        imag(SSva(ad.pq, pvq))  imag(SSvm(ad.pq, ad.pq))    ];
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
            ad = mm.get_userdata('aux_data');
            
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
            ad = mm.get_userdata('aux_data');
            Cp  = nm1.C([ad.pv; ad.pq], :);
            Cpp = nm2.C(ad.pq, :);
            Bp  = -imag( Cp  * Y1 * Cp' );
            Bpp = -imag( Cpp * Y2 * Cpp' );
            JJ = {Bp, Bpp};
        end

        function x = pf_gs_x_update(obj, x, f, mm, dm, mpopt);
            alg = mpopt.pf.alg;
            ad = mm.get_userdata('aux_data');

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
            ad = mm.get_userdata('aux_data');

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
                mm.userdata.aux_data = ad;
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
                    mm.userdata.aux_data = ad;
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
            mm.userdata.aux_data.vmpv = abs(v_(pv));

            x = [angle(v_(pvq)); abs(v_(pq))];
        end


        %%-----  CPF methods  -----
        function varargout = cpf_convert_x(obj, varargin)
            [varargout{1:nargout}] = obj.pf_convert_x(varargin{:});
        end

        function [f, J] = cpf_equations(obj, x, ad)
            %% index vector
            pvq = [ad.pv; ad.pq];

            b = [ real(ad.xfer(pvq)); imag(ad.xfer(ad.pq)) ];
            if nargout > 1
                [f, J] = obj.pf_node_balance_equations(x(1:end-1), ad);
                J = [J -b];
            else
                f = obj.pf_node_balance_equations(x(1:end-1), ad);
            end
            f = f - b * x(end);
        end

        function cpf_add_node_balance_constraints(obj, mm, dm, mpopt)
            ad = mm.get_userdata('aux_data');
            
            %% continuation power balance constraints
            fcn = @(x)cpf_equations(obj, x, ad);
            mm.add_nln_constraint({'Pmis', 'Qmis'}, [ad.npv+ad.npq;ad.npq], 1, fcn, []);
        end

        function ef = cpf_event_qlim(obj, cx, opt, mm, dm, mpopt)
            ad = mm.get_userdata('aux_data');

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
            ef = [v_Qmax; v_Qmin] * dm.baseMVA;
        end

        function [nx, cx, s] = cpf_callback_qlim(obj, k, nx, cx, px, s, opt, mm, dm, mpopt)
            %% skip if initialize, finalize or done
            if k <= 0 || s.done
                return;
            end

            %% handle event
            evnts = s.evnts;
            for i = 1:length(evnts)
                if strcmp(evnts(i).name, 'QLIM') && evnts(i).zero
                    efidx = evnts(i).idx;       %% event function index
                    [~, zi_min, zi_max] = obj.params_var('zi'); %% bounds on zi
                    if opt.verbose > 3
                        msg = sprintf('%s\n    ', evnts(i).msg);
                    else
                        msg = '';
                    end
                    for j = 1:length(efidx)
                        %% find index of z var
                        idx = efidx(j);         %% index of z var
                        if idx <= obj.nz
                            maxlim = 1;         %% Qmax violation
                            lim = zi_max(idx) * dm.baseMVA;
                            lim_type = 'Qmax';
                        else
                            idx = idx - obj.nz; %% correct index of z var
                            maxlim = 0;         %% Qmin violation
                            lim = zi_min(idx) * dm.baseMVA;
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
                        ad = mm.get_userdata('aux_data');
                        [v_, z_] = obj.cpf_convert_x(nx.x, ad);
                        z_(idx) = real(z_(idx)) + 1j * lim / dm.baseMVA;

                        %% change node type to PQ
                        obj.set_node_type_pq(dm, idx);

                        %% check for existence of remaining slack/PV bus
                        try
                            [ref, pv, pq] = obj.node_types(obj, dm);
                        catch
                            s.done = 1;
                            s.done_msg = 'No REF or PV nodes remaining.';
                        end

                        if ~s.done
                            %% get target case data, network models
                            pf = ad.target_pf;
                            dmt = pf.dm;
                            nmt = pf.nm;
                            
                            %% change node type in target case
                            nmt.set_node_type_pq(dmt, idx);

                            %% zero out Q transfer for bus
                            ad.xfer(idx) = 0;

                            %% if slack changed, zero out P transfer
                            if ref ~= ad.ref
                                ad.xfer(ref) = ad.xfer(ref) - real(ad.xfer(ref));
                            end

                            %% update aux data for math model
                            ad.ref  = ref;
                            ad.nref = length(ref);
                            ad.pv   = pv;
                            ad.npv  = length(pv);
                            ad.pq   = pq;
                            ad.npq  = length(pq);
                            mm.userdata.aux_data = ad;

                            nx.this_step = 0;
keyboard
                        %%    update bus types, in base and target cases
                        %%    select new slack bus if necessary
                        %%    (saving old 1st to be able to check for changes)
                        %%    update Q for newly limited inj in base & target (zero out transfer)
                        %%    if slack changed do the same for P at new slack bus
                        %%    nx.this_step = 0
                        end
                    end
                    s.evnts(i).msg = msg;
                end
            end


%             %% handle event
%             for i = 1:length(s.evnts)
%                 if strcmp(evnts(i).name, 'QLIM') && evnts(i).zero
% %                     %% get updated MPC, if necessary
% %                     if isempty(mpc)
% %                         d = cb_data;
% %                         if length(d.ref) ~= 1
% %                             error('cpf_qlim_event_cb: ''cpf.enforce_qlims'' option only valid for systems with exactly one REF bus')
% %                         end
% %                         mpc = cpf_current_mpc(d.mpc_base, d.mpc_target, ...
% %                             d.Ybus, d.Yf, d.Yt, d.ref, d.pv, d.pq, nx.V, nx.x(end), d.mpopt);
% %                         ng = size(mpc.gen, 1);
% %                         i2e_bus = mpc.order.bus.i2e;
% %                         i2e_gen = mpc.order.gen.i2e;
% %                     end
% 
%                     %% find the generator(s) and which lim(s)
%                     if opt.verbose > 3
%                         msg = sprintf('%s\n    ', evnts(i).msg);
%                     else
%                         msg = '';
%                     end
%                     ig = evnts(i).idx;
%                     for j = 1:length(ig)
%                         g = ig(j);                  %% index of gen of interest
%                         maxlim = 1;
%                         if g > obj.nz
%                             g = g - obj.nz;
%                             maxlim = 0;
%                         end
% %                         ib = mpc.gen(g, GEN_BUS);   %% corresponding bus index
%                         ib = g;
%                         [~, mn, mx] = obj.params_var('zi');
%                         if maxlim
%                             msg = sprintf('%sgen %d @ bus %d reached %g MVAr Qmax lim @ lambda = %.4g : bus %d converted to PQ', ...
%                                 msg, g, ib, mx(g), nx.x(end), ib);
%                         else
%                             msg = sprintf('%sgen %d @ bus %d reached %g MVAr Qmin lim @ lambda = %.4g : bus %d converted to PQ', ...
%                                 msg, g, ib, mn(g), nx.x(end), ib);
%                         end
% %                         if maxlim
% %                             msg = sprintf('%sgen %d @ bus %d reached %g MVAr Qmax lim @ lambda = %.4g : bus %d converted to PQ', ...
% %                                 msg, i2e_gen(g), i2e_bus(ib), mpc.gen(g, QMAX), nx.x(end), i2e_bus(ib));
% %                         else
% %                             msg = sprintf('%sgen %d @ bus %d reached %g MVAr Qmin lim @ lambda = %.4g : bus %d converted to PQ', ...
% %                                 msg, i2e_gen(g), i2e_bus(ib), mpc.gen(g, QMIN), nx.x(end), i2e_bus(ib));
% %                         end
% 
% %                         %% set Qg to exact limit and convert the generator's bus to PQ bus
% %                         if maxlim
% %                             mpc.gen(g, QG) = mpc.gen(g, QMAX);
% %                         else
% %                             mpc.gen(g, QG) = mpc.gen(g, QMIN);
% %                         end
% %                         mpc.bus(ib, BUS_TYPE) = PQ;
% % 
% %                         %% infeasibility check
% %                         on = find(mpc.gen(:, GEN_STATUS) > 0 & ...  %% which generators are on?
% %                                   mpc.bus(mpc.gen(:, GEN_BUS), BUS_TYPE) ~= PQ);  %% ... and are not PQ buses
% % 
% %                         if isempty(on)
% %                             s.done = 1;
% %                             s.done_msg = 'No REF or PV buses remaining.';
% %                         else
% %                             oldref = cb_data.ref;   %% save previous ref bus
% %                             [ref, pv, pq] = bustypes(mpc.bus, mpc.gen);
% %                             if oldref ~= ref        %% ref bus changed
% %                                 mpc.bus(ref, BUS_TYPE) = REF;
% %                             end
% % 
% %                             %% update new bus types, including in base and target cases
% %                             cb_data.ref = ref;
% %                             cb_data.pv  = pv;
% %                             cb_data.pq  = pq;
% %                             cb_data.mpc_base.bus(  :, BUS_TYPE) = mpc.bus(:, BUS_TYPE);
% %                             cb_data.mpc_target.bus(:, BUS_TYPE) = mpc.bus(:, BUS_TYPE);
% %             
% %                             %% update QG for Q limited gen in base and target
% %                             %% (no more reactive transfer for this gen)
% %                             cb_data.mpc_base.gen(  g, QG) = mpc.gen(g, QG);
% %                             cb_data.mpc_target.gen(g, QG) = mpc.gen(g, QG);
% % 
% %                             %% update PG for previous slack gen in base and target
% %                             %% (no more active transfer for this gen)
% %                             if oldref ~= ref
% %                                 cb_data.mpc_base.gen(  g, PG) = mpc.gen(g,PG);
% %                                 cb_data.mpc_target.gen(g, PG) = mpc.gen(g,PG);
% %                             end
% %                 
% %                             %% update functions
% %                             b = cb_data.mpc_base;
% %                             t = cb_data.mpc_target;
% %                             cb_data.Sbusb = @(Vm)makeSbus(b.baseMVA, b.bus, b.gen, d.mpopt, Vm);
% %                             cb_data.Sbust = @(Vm)makeSbus(t.baseMVA, t.bus, t.gen, d.mpopt, Vm);
%                 
%                             %% set size of next step to zero
%                             nx.this_step = 0;
% %                         end
%                     end
%                     evnts(i).msg = msg;
%                 end
%             end
        end

        %%-----  OPF methods  -----
        function opf_add_node_balance_constraints(obj, mm)
            %% power balance constraints
            nn = obj.node.N;            %% number of nodes
            fcn_mis = @(x)opf_power_balance_fcn(obj, obj.opf_convert_x(x));
            hess_mis = @(x, lam)opf_power_balance_hess(obj, ...
                obj.opf_convert_x(x), lam);
            mm.add_nln_constraint({'Pmis', 'Qmis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end

        function [lamP, lamQ] = opf_node_power_balance_prices(obj, mm)
            %% shadow prices on node power balance
            nne = mm.get_idx('nle');
            lambda = mm.soln.lambda;
            lamP = lambda.eqnonlin(nne.i1.Pmis:nne.iN.Pmis);
            lamQ = lambda.eqnonlin(nne.i1.Qmis:nne.iN.Qmis);
        end
    end     %% methods
end         %% classdef
