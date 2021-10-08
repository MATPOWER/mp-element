classdef mp_network_acci < mp_network_acc & mp_form_acci

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
        function [f, J] = pf_node_balance_equations(obj, x, ad)
            %% index vectors
            pvq = [ad.pv; ad.pq];
            pqv = [ad.pq; ad.pv];

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.pf_convert_x(x, ad, 1);

            %% incidence matrix
            C = obj.C;

            %% Jacobian
            if nargout > 1
                %% get port current injections with derivatives
                [I, dI.vr, dI.vi, dI.zr, dI.zi] = obj.port_inj_current([v_; z_], 1);
                dI.vr = C * dI.vr;
                dI.vi = C * dI.vi;
                dI.zr = C * dI.zr;
                dI.zi = C * dI.zi;
                %% derivatives of voltage magnitudes (for PV buses)
                nn = obj.node.N;
                dV2.vr = sparse(ad.pv, ad.pv, 2*real(v_(ad.pv)), nn, nn);
                dV2.vi = sparse(ad.pv, ad.pv, 2*imag(v_(ad.pv)), nn, nn);
                dV2.zr = sparse(nn, nn);
                dV2.zi = dV2.zr;
                JJ = cell(3, length(ad.var_map));

                for k = 1:length(ad.var_map)
                    m = ad.var_map{k};
                    name = m{1};
                    if ~isempty(m{2})       %% i1:iN
                        i1 = m{2};
                        iN = m{3};
                        JJ{1, k} = real(dI.(name)(pvq, i1:iN));
                        JJ{2, k} = imag(dI.(name)(pvq, i1:iN));
                        JJ{3, k} = dV2.(name)(ad.pv, i1:iN);
                    elseif isempty(m{4})    %% :
                        JJ{1, k} = real(dI.(name)(pvq, :));
                        JJ{2, k} = imag(dI.(name)(pvq, :));
                        JJ{3, k} = dV2.(name)(ad.pv, :);
                    else                    %% idx
                        idx = m{4};
                        JJ{1, k} = real(dI.(name)(pvq, idx));
                        JJ{2, k} = imag(dI.(name)(pvq, idx));
                        JJ{3, k} = dV2.(name)(ad.pv, idx);
                    end
                end
                J = vertcat( horzcat(JJ{1, :}), ...
                             horzcat(JJ{2, :}), ...
                             horzcat(JJ{3, :})  );
            else
                %% get port current injections (w/o derivatives)
                I = obj.port_inj_current([v_; z_], 1);
            end

            %% nodal power balance
            II = C * I;
            vmm = v_(ad.pv) .* conj(v_(ad.pv)) - ad.vr(ad.pv).^2 - ad.vi(ad.pv).^2;
            f = [real(II(pvq)); imag(II(pvq)); vmm];
        end

        %%-----  OPF methods  -----
        function [lam_p, lam_q] = opf_node_power_balance_prices(obj, mm)
            %% shadow prices on node power balance
            nne = mm.get_idx('nle');

            %% convert current balance shadow prices to equivalent lam_p and lam_q
            %% P + jQ = (Vr + jVi) * (M - jN)
            %% M = (Vr P + Vi Q) / (Vr^2 + Vi^2)
            %% N = (Vi P - Vr Q) / (Vr^2 + Vi^2)
            %% lam_p = df/dP = df/dM * dM/dP + df/dN + dN/dP
            %% lam_q = df/dQ = df/dM * dM/dQ + df/dN + dN/dQ
            V = obj.soln.v;
            lambda = mm.soln.lambda;
            VV = V ./ (V .* conj(V));   %% V / vm^2
            VVr = real(VV);
            VVi = imag(VV);
            lamM = lambda.eqnonlin(nne.i1.rImis:nne.iN.rImis);
            lamN = lambda.eqnonlin(nne.i1.iImis:nne.iN.iImis);
            lam_p = (VVr.*lamM + VVi.*lamN);
            lam_q = (VVi.*lamM - VVr.*lamN);
        end
    end     %% methods
end         %% classdef
