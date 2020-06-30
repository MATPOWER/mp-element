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
        function ad = power_flow_aux_data(obj, mpc, mpopt)
            %% call parent method
            ad = power_flow_aux_data@ac_aggregate(obj, mpc, mpopt);

            %% build additional aux data
            g = obj.mpe_by_name('gen');
            i = obj.mpe_z_map(obj.mpe_map.gen, :);  %% 1st-last z-idx for gens
            N = g.C(ad.pv, :) * g.N;%% z coefficients for all gens @ PV nodes
            [ii, jj, ss] = find(N); %% deconstruct and recreate with
            [~, ia] = unique(ii);   %% only 1st non-zero in each row
            N = sparse(ii(ia), jj(ia), ss(ia), ad.npv, size(N, 2));
            j = find(any(N, 1));    %% indices of PV node gen z-vars (in gen z)
            k = j + i(1) - 1;       %% indices of PV node gen z-vars (in sys z)
            N = N(:,j) * g.D(k,j)'; %% coefficients for zi(k)

            %% save additional aux data
            ad.N = N;               %% coefficients for zi(k) corrsp to PV node gens
            ad.invN = inv(N);       %% inverse of N, typically equal to N = -eye()
            ad.k = k;               %% indices of PV node gen z-vars (in sys z)
        end

        function add_pf_vars(obj, asm, om, ad, mpc, mpopt)
            %% get model variables
            vvars = obj.model_vvars();

            %% index vectors
            pqv = [ad.pq; ad.pv];

            %% reactive injections
            v_ = ad.v1 + 1j * ad.v2;
            z_ = ad.zr + 1j * ad.zi;
            Qpv = obj.C(ad.pv, :) * imag( obj.port_inj_power([v_; z_], 1) );
            Qg_pv = Qpv - ad.N * ad.zi(ad.k);
            om.add_var('Qg_pv', ad.npv, Qg_pv);

            %% voltage real part
            st = obj.(vvars{1});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    om.add_var(name, ad.npq+ad.npv, d.v0.(name)(pqv), d.vl.(name)(pqv), d.vu.(name)(pqv));
                else
                    error('handling of indexed sets not implmented here (yet)');
                end
            end

            %% voltage imaginary part
            st = obj.(vvars{2});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    om.add_var(name, ad.npq+ad.npv, d.v0.(name)(pqv), d.vl.(name)(pqv), d.vu.(name)(pqv));
                else
                    error('handling of indexed sets not implmented here (yet)');
                end
            end
        end

        function [v_, z_] = pfx2vz(obj, x, ad)
            %% update v_, z_ from x
            pqv = [ad.pq; ad.pv];
            iN = ad.npv;                            Qg_pv   = x(1:iN);
            i1 = iN+1;  iN = iN + ad.npv + ad.npq;  ad.v1(pqv) = x(i1:iN);  %% vr
            i1 = iN+1;  iN = iN + ad.npv + ad.npq;  ad.v2(pqv) = x(i1:iN);  %% vi
            v_ = ad.v1 + 1j * ad.v2;
            ad.zi(ad.k) = -ad.N \ Qg_pv;
            z_ = ad.zr + 1j * ad.zi;
        end

        function [F, J] = power_flow_equations(obj, x, ad)
            %% index vectors
            pvq = [ad.pv; ad.pq];
            pqv = [ad.pq; ad.pv];

            %% update model state ([v_; z_]) from power flow state (x)
            [v_, z_] = obj.pfx2vz(x, ad);

            %% incidence matrix
            C = obj.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                [I, Iva, Ivm, Izr, Izi] = obj.port_inj_current([v_; z_], 1);
                dV2_dVr = sparse(1:ad.npv, ad.npq+(1:ad.npv), 2*real(v_(ad.pv)), ad.npv, ad.npv+ad.npq);
                dV2_dVi = sparse(1:ad.npv, ad.npq+(1:ad.npv), 2*imag(v_(ad.pv)), ad.npv, ad.npv+ad.npq);

                IIva = C * Iva;
                IIvm = C * Ivm;
                IIzi = C * Izi;
                dImis_dQg = -IIzi(:, ad.k) * ad.invN;

                J = [   real(dImis_dQg(pvq, :)) real(IIva(pvq, pqv)) real(IIvm(pvq, pqv));
                        imag(dImis_dQg(pvq, :)) imag(IIva(pvq, pqv)) imag(IIvm(pvq, pqv));
                        sparse(ad.npv, ad.npv) dV2_dVr dV2_dVi ];
            else
                %% get port power injections (w/o derivatives)
                I = obj.port_inj_current([v_; z_], 1);
            end

            %% nodal power balance
            II = C * I;
            vmm = v_(ad.pv) .* conj(v_(ad.pv)) - ad.v1(ad.pv).^2 - ad.v2(ad.pv).^2;
            F = [real(II(pvq)); imag(II(pvq)); vmm];
        end

        function add_pf_node_balance_constraints(obj, om, ad)
            %% power balance constraints
            npvq = ad.npv+ad.npq;
            fcn = @(x)power_flow_equations(obj, x, ad);
            om.add_nln_constraint({'Irmis', 'Iimis', 'Vmis'}, [npvq;npvq;ad.npv], 1, fcn, []);
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
