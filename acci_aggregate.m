classdef acci_aggregate < acc_aggregate% & acci_model

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
        function ad = power_flow_aux_data(obj, vr, vi, zr, zi, t)
            g = obj.mpe_by_name('gen');
            i = obj.mpe_z_map(obj.mpe_map.gen, :);  %% 1st-last z-idx for gens
            N = g.C(t.pv, :) * g.N; %% z coefficients for all gens @ PV nodes
            [ii, jj, ss] = find(N); %% deconstruct and recreate with
            [~, ia] = unique(ii);   %% only 1st non-zero in each row
            N = sparse(ii(ia), jj(ia), ss(ia), t.npv, size(N, 2));
            j = find(any(N, 1));    %% indices of PV node gen z-vars (in gen z)
            k = j + i(1) - 1;       %% indices of PV node gen z-vars (in sys z)
            N = N(:,j) * g.D(k,j)'; %% coefficients for zi(k)
            ad  = struct( ...
                'N', N, ...         %% coefficients for zi(k) corrsp to PV node gens
                'invN', inv(N), ... %% inverse of N, typically equal to N = -eye()
                'k', k ...          %% indices of PV node gen z-vars (in sys z)
            );
        end

        function x = vz2pfx(obj, vr, vi, zr, zi, t, ad)
            %% update x from vr, vi, zr, zi
            pqv = [t.pq; t.pv];
            v_ = vr + 1j * vi;
            z_ = zr + 1j * zi;
            Qpv = obj.C(t.pv, :) * imag( obj.port_inj_power([v_; z_], 1) );
            Qg_pv = Qpv - ad.N * zi(ad.k);
            x = [Qg_pv; vr(pqv); vi(pqv)];
        end

        function [v_, z_] = pfx2vz(obj, x, vr, vi, zr, zi, t, ad)
            %% update v_, z_ from x
            pqv = [t.pq; t.pv];
            iN = t.npv;                             Qg_pv   = x(1:iN);
            i1 = iN+1;  iN = iN + t.npv + t.npq;    vr(pqv) = x(i1:iN);
            i1 = iN+1;  iN = iN + t.npv + t.npq;    vi(pqv) = x(i1:iN);
            v_ = vr + 1j * vi;
            zi(ad.k) = -ad.N \ Qg_pv;
            z_ = zr + 1j * zi;
        end

        function [F, J] = power_flow_equations(obj, x, vr, vi, zr, zi, t, ad)
            %% index vectors
            pvq = [t.pv; t.pq];
            pqv = [t.pq; t.pv];

            %% update model state ([v_; z_]) from power flow state (x)
            [v_, z_] = pfx2vz(obj, x, vr, vi, zr, zi, t, ad);

            %% incidence matrix
            C = obj.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                [I, Iva, Ivm, Izr, Izi] = obj.port_inj_current([v_; z_], 1);
                dV2_dVr = sparse(1:t.npv, t.npq+(1:t.npv), 2*real(v_(t.pv)), t.npv, t.npv+t.npq);
                dV2_dVi = sparse(1:t.npv, t.npq+(1:t.npv), 2*imag(v_(t.pv)), t.npv, t.npv+t.npq);

                IIva = C * Iva;
                IIvm = C * Ivm;
                IIzi = C * Izi;
                dImis_dQg = -IIzi(:, ad.k) * ad.invN;

                J = [   real(dImis_dQg(pvq, :)) real(IIva(pvq, pqv)) real(IIvm(pvq, pqv));
                        imag(dImis_dQg(pvq, :)) imag(IIva(pvq, pqv)) imag(IIvm(pvq, pqv));
                        sparse(t.npv, t.npv) dV2_dVr dV2_dVi ];
            else
                %% get port power injections (w/o derivatives)
                I = obj.port_inj_current([v_; z_], 1);
            end

            %% nodal power balance
            II = C * I;
            vmm = v_(t.pv) .* conj(v_(t.pv)) - vr(t.pv).^2 - vi(t.pv).^2;
            F = [real(II(pvq)); imag(II(pvq)); vmm];
        end


        %%-----  OPF methods  -----
        function add_opf_node_balance_constraints(obj, om)
            %% power balance constraints
            nn = obj.node.N;            %% number of nodes
            fcn_mis = @(x)opf_current_balance_fcn(obj, x);
            hess_mis = @(x, lam)opf_current_balance_hess(obj, x, lam);
            om.add_nln_constraint({'rImis', 'iImis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end
    end     %% methods
end         %% classdef
