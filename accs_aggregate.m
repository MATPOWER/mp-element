classdef accs_aggregate < acc_aggregate% & accs_model

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
        function add_opf_node_balance_constraints(obj, om)
            %% power balance constraints
            nn = obj.node.N;            %% number of nodes
            fcn_mis = @(x)opf_power_balance_fcn(obj, x);
            hess_mis = @(x, lam)opf_power_balance_hess(obj, x, lam);
            om.add_nln_constraint({'Pmis', 'Qmis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end


        %%-----  PF methods  -----
        function x = vz2pfx(obj, vr, vi, zr, zi, t, ad)
            %% update x from vr, vi, zr, zi
            pqv = [t.pq; t.pv];
            x = [vr(pqv); vi(pqv)];
        end

        function [v_, z_] = pfx2vz(obj, x, vr, vi, zr, zi, t, ad)
            %% update v_, z_ from x
            pqv = [t.pq; t.pv];
            vr(pqv) = x(1:t.npv+t.npq);
            vi(pqv) = x(t.npv+t.npq+1:end);
            v_ = vr + 1j * vi;
            z_ = zr + 1j * zi;
        end

        function [F, J] = power_flow_equations(obj, x, vr, vi, zr, zi, t, ad)
            %% index vector
            pqv = [t.pq; t.pv];

            %% update model state ([v_; z_]) from power flow state (x)
            [v_, z_] = pfx2vz(obj, x, vr, vi, zr, zi, t);

            %% incidence matrix
            C = obj.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                [S, Svr, Svi] = obj.port_inj_power([v_; z_], 1);
                dV2_dVr = sparse(1:t.npv, t.npq+(1:t.npv), 2*real(v_(t.pv)), t.npv, t.npv+t.npq);
                dV2_dVi = sparse(1:t.npv, t.npq+(1:t.npv), 2*imag(v_(t.pv)), t.npv, t.npv+t.npq);

                SSvr = C * Svr;
                SSvi = C * Svi;
                J = [   real(SSvr(pqv,  pqv)) real(SSvi(pqv,  pqv));
                        imag(SSvr(t.pq, pqv)) imag(SSvi(t.pq, pqv));
                        dV2_dVr dV2_dVi ];
            %     SSzr = C * Szr;
            %     SSzi = C * Szi;
            %     J = [
            %         real(SSvr(pvq,  pvq)) real(SSvi(pvq,  t.pq)) real(SSzr(pvq,  :)) real(SSzi(pvq,  :));
            %         imag(SSvr(t.pq, pvq)) imag(SSvi(t.pq, t.pq)) imag(SSzr(t.pq, :)) imag(SSzi(t.pq, :))  ];
            else
                %% get port power injections (w/o derivatives)
                S = obj.port_inj_power([v_; z_], 1);
            end

            %% nodal power balance + voltage magnitude mismatch
            SS = C * S;
            vmm = v_(t.pv) .* conj(v_(t.pv)) - vr(t.pv).^2 - vi(t.pv).^2;
            F = [real(SS(pqv)); imag(SS(t.pq)); vmm];
        end
    end     %% methods
end         %% classdef
