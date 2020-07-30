classdef acps_aggregate < acp_aggregate% & acps_model

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
        function add_pf_vars(obj, asm, om, mpc, mpopt)
            %% get model variables
            vvars = obj.model_vvars();

            %% index vectors
            ad = om.get_userdata('power_flow_aux_data');
            pvq = [ad.pv; ad.pq];

            %% voltage angles
            st = obj.(vvars{1});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    om.add_var(name, ad.npv+ad.npq, d.v0.(name)(pvq), d.vl.(name)(pvq), d.vu.(name)(pvq));
                else
                    error('handling of indexed sets not implmented here (yet)');
                end
            end

            %% voltage magnitudes
            st = obj.(vvars{2});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    om.add_var(name, ad.npq, d.v0.(name)(ad.pq), d.vl.(name)(ad.pq), d.vu.(name)(ad.pq));
                else
                    error('handling of indexed sets not implmented here (yet)');
                end
            end
        end

        function [v_, z_] = pfx2vz(obj, x, ad)
            %% update v_, z_ from x
            ad.v1([ad.pv; ad.pq]) = x(1:ad.npv+ad.npq);         %% va
            ad.v2(ad.pq)          = x(ad.npv+ad.npq+1:end);     %% vm
            v_ = ad.v2 .* exp(1j * ad.v1);
            z_ = ad.zr + 1j * ad.zi;
        end

        function [f, J] = power_flow_equations(obj, x, ad, fdpf)
            %% index vector
            pvq = [ad.pv; ad.pq];

            %% update model state ([v_; z_]) from power flow state (x)
            [v_, z_] = obj.pfx2vz(x, ad);

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

        function add_pf_node_balance_constraints(obj, om, mpc, mpopt)
            alg = mpopt.pf.alg;
            ad = om.get_userdata('power_flow_aux_data');
            
            %% power balance constraints
            switch alg
                case  {'FDXB', 'FDBX'}
                    fcn = @(x)power_flow_equations(obj, x, ad, 1);
                otherwise
                    fcn = @(x)power_flow_equations(obj, x, ad);
            end
            om.add_nln_constraint({'Pmis', 'Qmis'}, [ad.npv+ad.npq;ad.npq], 1, fcn, []);
        end

        function JJ = fd_jac_approx(obj, om, mpc, mpopt)
            alg = mpopt.pf.alg;
            
            %% define named indices into bus, branch matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

            %% modify data model to form Bp (B prime)
            mpc1 = mpc;
            mpc1.bus(:, BS) = 0;            %% zero out shunts at buses
            mpc2 = mpc1;
            mpc1.branch(:, BR_B) = 0;       %% zero out line charging shunts
            mpc1.branch(:, TAP) = 1;        %% cancel out taps
            if strcmp(alg, 'FDXB')          %% if XB method
                mpc1.branch(:, BR_R) = 0;   %% zero out line resistance
            end

            %% modify data model to form Bpp (B double prime)
            mpc2.branch(:, SHIFT) = 0;      %% zero out phase shifters
            if strcmp(alg, 'FDBX')          %% if BX method
                mpc2.branch(:, BR_R) = 0;   %% zero out line resistance
            end

            %% build network models and get admittance matrices
            nm1 = feval(class(obj)).create_model(mpc1, mpopt);
            nm2 = feval(class(obj)).create_model(mpc2, mpopt);
            [Y1, L, M] = nm1.get_params([], {'Y', 'L', 'M'});
            Y2 = nm2.get_params();
            if any(any(L)) || any(any(M))
                error('acps_fdpf_aggregate/df_jac_approx: fast-decoupled Jacobian approximation not implemented for models with non-zero L and/or M matrices.')
            end

            %% form full Bp and Bpp matrices
            ad = om.get_userdata('power_flow_aux_data');
            Cp  = nm1.C([ad.pv; ad.pq], :);
            Cpp = nm2.C(ad.pq, :);
            Bp  = -imag( Cp  * Y1 * Cp' );
            Bpp = -imag( Cpp * Y2 * Cpp' );
            JJ = {Bp, Bpp};
        end

        function x = gs_x_update(obj, x, f, om, mpc, mpopt);
            alg = mpopt.pf.alg;
            ad = om.get_userdata('power_flow_aux_data');

            %% get model state ([v_; z_]) from power flow state (x)
            [v_, z_] = obj.pfx2vz(x, ad);

            C = obj.C;
            Y = C * obj.get_params([], 'Y') * C';
            SS = zeros(size(C, 1), 1);
            SS(ad.pv) = f(1:ad.npv);
            SS(ad.pq) = f(ad.npv+1:ad.npv+ad.npq) + ...
                    1j * f(ad.npv+ad.npq+1:ad.npv+2*ad.npq);
%            SS = C * obj.port_inj_power([v_; z_], 1);
            S0 = v_ .* conj(Y * v_) - SS;

            %% update voltage
            %% at PQ buses
            for k = ad.pq'
                v_(k) = v_(k) + (conj(S0(k)/v_(k)) - Y(k,:) * v_) / Y(k, k);
            end

            %% at PV buses
            if ad.npv
                for k = ad.pv'
                    S0(k) = real(S0(k)) + 1j * imag( v_(k) * conj(Y(k,:) * v_) );
                    v_(k) = v_(k) + (conj(S0(k)/v_(k)) - Y(k,:) * v_) / Y(k, k);
                end
            end

            x = [angle(v_([ad.pv; ad.pq])); abs(v_(ad.pq))];
        end


        %%-----  OPF methods  -----
        function add_opf_node_balance_constraints(obj, om)
            %% power balance constraints
            nn = obj.node.N;            %% number of nodes
            fcn_mis = @(x)opf_power_balance_fcn(obj, x);
            hess_mis = @(x, lam)opf_power_balance_hess(obj, x, lam);
            om.add_nln_constraint({'Pmis', 'Qmis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end
    end     %% methods
end         %% classdef
