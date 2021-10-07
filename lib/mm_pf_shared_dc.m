classdef mm_pf_shared_dc < mm_pf_shared

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end
    
    methods
        function ad = pf_aux_data(obj, nm, dm, mpopt)
            %% call parent
            ad = pf_aux_data@mm_pf_shared(obj, nm, dm, mpopt);

            %% get parameters
            [B, K, p] = nm.get_params();

            %% add DC model parameters
            ad.B = nm.C * B * nm.C';
            ad.Pbus = -(nm.C * K * nm.D' * ad.z + nm.C * p);
        end

        function obj = add_pf_system_vars(obj, nm, dm, mpopt)
            %% get model variables
            vvars = nm.model_vvars();

            %% index vectors
            ad = obj.aux_data;
            pvq = [ad.pv; ad.pq];

            %% voltage angles
            st = nm.(vvars{1});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    obj.add_var(name, ad.npv+ad.npq, d.v0.(name)(pvq), d.vl.(name)(pvq), d.vu.(name)(pvq));
                else
                    error('mp_math_pf/add_pf_system_vars: handling of indexed sets not implmented here (yet)');
                end
            end
        end
    end     %% methods
end         %% classdef
