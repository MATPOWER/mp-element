classdef (Abstract) mm_shared_pfcpf_accs < mp.mm_shared_pfcpf_acc

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        function obj = add_system_vars_pf(obj, nm, dm, mpopt)
            %% get model variables
            vvars = nm.model_vvars();

            %% voltage real part
            obj.add_system_varset_pf(nm, vvars{1}, 'pq');
            obj.add_system_varset_pf(nm, vvars{1}, 'pv');

            %% voltage imaginary part
            obj.add_system_varset_pf(nm, vvars{2}, 'pq');
            obj.add_system_varset_pf(nm, vvars{2}, 'pv');
        end

        function [f, J] = node_balance_equations(obj, x, nm)
            %% index vector
            ad = obj.aux_data;
            pqv = [ad.pq; ad.pv];

            %% update network model state ([v_; z_]) from math model state (x)
            [v_, z_] = obj.convert_x_m2n(x, nm, 1);

            %% incidence matrix
            C = nm.C;

            %% Jacobian
            if nargout > 1
                %% get port power injections with derivatives
                var_names = cellfun(@(x)x{1}, ad.var_map, 'UniformOutput', false);
                dz = any(strcmp(var_names, 'zr')) || ...
                     any(strcmp(var_names, 'zi'));
                if dz
                    [S, dS.vr, dS.vi, dS.zr, dS.zi] = nm.port_inj_power([v_; z_], 1);
                else
                    [S, dS.vr, dS.vi] = nm.port_inj_power([v_; z_], 1);
                end
                dS.vr = C * dS.vr;
                dS.vi = C * dS.vi;
                if dz
                    dS.zr = C * dS.zr;
                    dS.zi = C * dS.zi;
                end

                %% derivatives of voltage magnitudes (for PV buses)
                nn = nm.node.N;
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
                        JJ{1, k} = real(dS.(name)(pqv,   i1:iN));
                        JJ{2, k} = imag(dS.(name)(ad.pq, i1:iN));
                        JJ{3, k} = dV2.(name)(ad.pv, i1:iN);
                    elseif isempty(m{4})    %% :
                        JJ{1, k} = real(dS.(name)(pqv,   :));
                        JJ{2, k} = imag(dS.(name)(ad.pq, :));
                        JJ{3, k} = dV2.(name)(ad.pv, :);
                    else                    %% idx
                        idx = m{4};
                        JJ{1, k} = real(dS.(name)(pqv,   idx));
                        JJ{2, k} = imag(dS.(name)(ad.pq, idx));
                        JJ{3, k} = dV2.(name)(ad.pv, idx);
                    end
                end
                J = vertcat( horzcat(JJ{1, :}), ...
                             horzcat(JJ{2, :}), ...
                             horzcat(JJ{3, :})  );
            else
                %% get port power injections (w/o derivatives)
                S = nm.port_inj_power([v_; z_], 1);
            end

            %% nodal power balance + voltage magnitude mismatch
            SS = C * S;
            vmm = v_(ad.pv) .* conj(v_(ad.pv)) - ad.vr(ad.pv).^2 - ad.vi(ad.pv).^2;
            f = [real(SS(pqv)); imag(SS(ad.pq)); vmm];
        end
    end     %% methods
end         %% classdef
