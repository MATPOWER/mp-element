classdef mp_network_acpi < mp_network_acp% & mp_form_acpi

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
        function ad = power_flow_aux_data(obj, dm, mpopt)
            %% call parent method
            ad = power_flow_aux_data@mp_network_ac(obj, dm, mpopt);

            %% build additional aux data
            g = obj.elm_by_name('gen');
            i = obj.nme_z_map(obj.elm_map.gen, :);  %% 1st-last z-idx for gens
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

        function add_pf_vars(obj, nm, om, dm, mpopt)
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
                    error('mp_network_acpi/add_pf_vars: handling of indexed sets not implmented here (yet)');
                end
            end

            %% reactive injections
            v_ = ad.v1 + 1j * ad.v2;
            z_ = ad.zr + 1j * ad.zi;
            Qpv = obj.C(ad.pv, :) * imag( obj.port_inj_power([v_; z_], 1) );
            Qg_pv = Qpv - ad.N * ad.zi(ad.k);
            om.add_var('Qg_pv', ad.npv, Qg_pv);

            %% voltage magnitudes
            st = obj.(vvars{2});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    om.add_var(name, ad.npq, d.v0.(name)(ad.pq), d.vl.(name)(ad.pq), d.vu.(name)(ad.pq));
                else
                    error('mp_network_acpi/add_pf_vars: handling of indexed sets not implmented here (yet)');
                end
            end
        end

        function [v_, z_] = pf_convert_x(obj, x, ad)
            %% update v_, z_ from x
            iN = ad.npv + ad.npq;           ad.v1([ad.pv; ad.pq]) = x(1:iN);%% va
            i1 = iN+1;  iN = iN + ad.npv;   Qg_pv = x(i1:iN);
            i1 = iN+1;  iN = iN + ad.npq;   ad.v2(ad.pq) = x(i1:iN);        %% vm
            v_ = ad.v2 .* exp(1j * ad.v1);
            ad.zi(ad.k) = -ad.N \ Qg_pv;
            z_ = ad.zr + 1j * ad.zi;
        end

        function [f, J] = power_flow_equations(obj, x, ad)
            %% index vector
            pvq = [ad.pv; ad.pq];

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.pf_convert_x(x, ad);

            %% incidence matrix
            C = obj.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                [I, Iva, Ivm, Izr, Izi] = obj.port_inj_current([v_; z_], 1);

                IIva = C * Iva;
                IIvm = C * Ivm;
                IIzi = C * Izi;
                IIvm(:, ad.pv) = -IIzi(:, ad.k) * ad.invN;  %% dImis_dQg

                J = [   real(IIva(pvq, pvq))    real(IIvm(pvq, pvq));
                        imag(IIva(pvq, pvq))    imag(IIvm(pvq, pvq))    ];
            else
                %% get port power injections (w/o derivatives)
                I = obj.port_inj_current([v_; z_], 1);
            end

            %% nodal power balance
            II = C * I;
            f = [real(II(pvq)); imag(II(pvq))];
        end

        function add_pf_node_balance_constraints(obj, om, dm, mpopt)
            %% power balance constraints
            ad = om.get_userdata('power_flow_aux_data');
            npvq = ad.npv+ad.npq;
            fcn = @(x)power_flow_equations(obj, x, ad);
            om.add_nln_constraint({'Irmis', 'Iimis'}, [npvq;npvq], 1, fcn, []);
        end


        %%-----  OPF methods  -----
        function add_opf_node_balance_constraints(obj, om)
            %% power balance constraints
            nn = obj.node.N;            %% number of nodes
            fcn_mis = @(x)opf_current_balance_fcn(obj, obj.opf_convert_x(x));
            hess_mis = @(x, lam)opf_current_balance_hess(obj, ...
                obj.opf_convert_x(x), lam);
            om.add_nln_constraint({'rImis', 'iImis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end
    end     %% methods
end         %% classdef
