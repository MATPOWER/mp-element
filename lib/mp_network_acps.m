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

        function opt = cpf_solve_opts_warmstart(obj, opt, mm)
            ws = opt.warmstart;
            ad = mm.get_userdata('aux_data');

            %% set starting point to where previous continuation left off
            opt.x0 = [angle(ws.V([ad.pv; ad.pq])); abs(ws.V(ad.pq)); ws.lam];
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
                    ad = mm.get_userdata('aux_data');
                    if ad.nref ~= 1
                        error('mp_network_acps/cpf_callback_qlim: ''cpf.enforce_qlims'' option only valid for systems with exactly one REF bus');
                    end

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
                        [v_, z_] = obj.cpf_convert_x(nx.x, ad);
                        z_(idx) = real(z_(idx)) + 1j * lim / dm.baseMVA;

                        %% change node type to PQ
                        obj.set_node_type_pq(dm, nidx);

                        %% check for existence of remaining slack/PV bus
                        try
                            %% potentially pick new reference bus
                            [ref, pv, pq] = obj.node_types(obj, dm);
                        catch
                            s.done = 1;
                            s.done_msg = 'No REF or PV nodes remaining.';
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
                        s.warmstart = struct('nmt', nmt, 'dmt', dmt);
                    end
                    s.evnts(i).msg = msg;
                end
            end
        end
                    s.evnts(i).msg = msg;
                end
            end
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
