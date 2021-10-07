classdef mm_pf_shared_ac_i < mm_pf_shared

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

            %% build additional aux data
            N = nm.C(ad.pv, :) * nm.N;      %% z coefficients for z @ PV nodes
            [ii, jj, ~] = find(N == -1);    %% find element representing 1st
            [~, ia] = unique(ii, 'first');  %% direct injection at each PV node
            [~, ib] = sort(ii(ia));         %% sort by PV node

            %% save additional aux data
            ad.zi_idx = jj(ia(ib));         %% z-var index corresp to PV nodes
        end
    end     %% methods
end         %% classdef
