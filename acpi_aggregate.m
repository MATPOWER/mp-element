classdef acpi_aggregate < acp_aggregate% & acpi_model

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
        function ad = power_flow_aux_data(obj, va, vm, zr, zi, t)
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

        function x = vz2pfx(obj, va, vm, zr, zi, t, ad)
            %% update x from va, vm, zr, zi
            v_ = vm .* exp(1j * va);
            z_ = zr + 1j * zi;
            Qpv = obj.C(t.pv, :) * imag( obj.port_inj_power([v_; z_], 1) );
            Qg_pv = Qpv - ad.N * zi(ad.k);
            x = [va([t.pv; t.pq]); Qg_pv; vm(t.pq)];
        end

        function [v_, z_] = pfx2vz(obj, x, va, vm, zr, zi, t, ad)
            %% update v_, z_ from x
            iN = t.npv + t.npq;             va([t.pv; t.pq]) = x(1:iN);
            i1 = iN+1;  iN = iN + t.npv;    Qg_pv = x(i1:iN);
            i1 = iN+1;  iN = iN + t.npq;    vm(t.pq) = x(i1:iN);
            v_ = vm .* exp(1j * va);
            zi(ad.k) = -ad.N \ Qg_pv;
            z_ = zr + 1j * zi;
        end

        function [F, J] = power_flow_equations(obj, x, va, vm, zr, zi, t, ad)
            %% index vector
            pvq = [t.pv; t.pq];

            %% update model state ([v_; z_]) from power flow state (x)
            [v_, z_] = pfx2vz(obj, x, va, vm, zr, zi, t, ad);

            %% incidence matrix
            C = obj.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                [I, Iva, Ivm, Izr, Izi] = obj.port_inj_current([v_; z_], 1);

                IIva = C * Iva;
                IIvm = C * Ivm;
                IIzi = C * Izi;
                IIvm(:, t.pv) = -IIzi(:, ad.k) * ad.invN;   %% dImis_dQg

                J = [   real(IIva(pvq, pvq)) real(IIvm(pvq, pvq));
                        imag(IIva(pvq, pvq)) imag(IIvm(pvq, pvq))  ];
            else
                %% get port power injections (w/o derivatives)
                I = obj.port_inj_current([v_; z_], 1);
            end

            %% nodal power balance
            II = C * I;
            F = [real(II(pvq)); imag(II(pvq))];
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
