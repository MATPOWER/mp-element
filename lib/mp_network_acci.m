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
