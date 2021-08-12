classdef mp_network_accs < mp_network_acc & mp_form_accs

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
        function obj = pf_add_vars(obj, mm, nm, dm, mpopt)
            %% get model variables
            vvars = obj.model_vvars();

            %% index vectors
            ad = mm.get_userdata('aux_data');
            pqv = [ad.pq; ad.pv];

            %% voltage real part
            st = obj.(vvars{1});
            for k = 1:st.NS
                name = st.order(k).name;
                pq = ad.node_type_by_elm(k).pq;
                npq = length(pq);
                if isempty(st.order(k).idx)
                    d = st.data;
                    mm.add_var([name '_pq'], npq, d.v0.(name)(pq), d.vl.(name)(pq), d.vu.(name)(pq));
                else
                    error('mp_network_accs/pf_add_vars: handling of indexed sets not implmented here (yet)');
                end
            end
            for k = 1:st.NS
                name = st.order(k).name;
                pv = ad.node_type_by_elm(k).pv;
                npv = length(pv);
                if isempty(st.order(k).idx)
                    d = st.data;
                    mm.add_var([name '_pv'], npv, d.v0.(name)(pv), d.vl.(name)(pv), d.vu.(name)(pv));
                else
                    error('mp_network_accs/pf_add_vars: handling of indexed sets not implmented here (yet)');
                end
            end

            %% voltage imaginary part
            st = obj.(vvars{2});
            for k = 1:st.NS
                name = st.order(k).name;
                pq = ad.node_type_by_elm(k).pq;
                npq = length(pq);
                if isempty(st.order(k).idx)
                    d = st.data;
                    mm.add_var([name '_pq'], npq, d.v0.(name)(pq), d.vl.(name)(pq), d.vu.(name)(pq));
                else
                    error('mp_network_accs/pf_add_vars: handling of indexed sets not implmented here (yet)');
                end
            end
            for k = 1:st.NS
                name = st.order(k).name;
                pv = ad.node_type_by_elm(k).pv;
                npv = length(pv);
                if isempty(st.order(k).idx)
                    d = st.data;
                    mm.add_var([name '_pv'], npv, d.v0.(name)(pv), d.vl.(name)(pv), d.vu.(name)(pv));
                else
                    error('mp_network_accs/pf_add_vars: handling of indexed sets not implmented here (yet)');
                end
            end
        end

        function [vx_, z_, x_] = pf_convert_x(obj, mmx, ad, only_v)
            %% x = obj.pf_convert(mmx, ad)
            %% [v, z] = obj.pf_convert(mmx, ad)
            %% [v, z, x] = obj.pf_convert(mmx, ad)
            %% ... = obj.pf_convert(mmx, ad, only_v)

            %% update v_, z_ from mmx
            pqv = [ad.pq; ad.pv];
            ad.v1(pqv) = mmx(1:ad.npv+ad.npq);                      %% vr
            ad.v2(pqv) = mmx(ad.npv+ad.npq+1:2*ad.npv+2*ad.npq);    %% vi
            vx_ = ad.v1 + 1j * ad.v2;
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

        function [f, J] = pf_node_balance_equations(obj, x, ad)
            %% index vector
            pqv = [ad.pq; ad.pv];

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.pf_convert_x(x, ad, 1);

            %% incidence matrix
            C = obj.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                [S, Svr, Svi] = obj.port_inj_power([v_; z_], 1);
                dV2_dVr = sparse(1:ad.npv, ad.npq+(1:ad.npv), 2*real(v_(ad.pv)), ad.npv, ad.npv+ad.npq);
                dV2_dVi = sparse(1:ad.npv, ad.npq+(1:ad.npv), 2*imag(v_(ad.pv)), ad.npv, ad.npv+ad.npq);

                SSvr = C * Svr;
                SSvi = C * Svi;
                J = [   real(SSvr(pqv,  pqv))   real(SSvi(pqv,  pqv));
                        imag(SSvr(ad.pq, pqv))  imag(SSvi(ad.pq, pqv));
                        dV2_dVr dV2_dVi ];
            else
                %% get port power injections (w/o derivatives)
                S = obj.port_inj_power([v_; z_], 1);
            end

            %% nodal power balance + voltage magnitude mismatch
            SS = C * S;
            vmm = v_(ad.pv) .* conj(v_(ad.pv)) - ad.v1(ad.pv).^2 - ad.v2(ad.pv).^2;
            f = [real(SS(pqv)); imag(SS(ad.pq)); vmm];
        end

        function obj = pf_add_node_balance_constraints(obj, mm, dm, mpopt)
            %% power balance constraints
            ad = mm.get_userdata('aux_data');
            fcn = @(x)pf_node_balance_equations(obj, x, ad);
            mm.add_nln_constraint({'Pmis', 'Qmis', 'Vmis'}, [ad.npv+ad.npq;ad.npq;ad.npv], 1, fcn, []);
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

        function [lam_p, lam_q] = opf_node_power_balance_prices(obj, mm)
            %% shadow prices on node power balance
            nne = mm.get_idx('nle');
            lambda = mm.soln.lambda;
            lam_p = lambda.eqnonlin(nne.i1.Pmis:nne.iN.Pmis);
            lam_q = lambda.eqnonlin(nne.i1.Qmis:nne.iN.Qmis);
        end
    end     %% methods
end         %% classdef
