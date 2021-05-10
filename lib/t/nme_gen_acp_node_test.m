classdef nme_gen_acp_node_test < nme_gen_acp

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end
%     
    methods
        function idx = node_indices(obj, nm, dm, dme)
            bus_dme = dm.elm_by_name('bus_ld');
            nidx = nm.get_node_idx('bus_ld');   %% node indices for 'bus'
            b = dme.bus(dme.on);
            bidx = bus_dme.i2on(b);         %% online bus indices for gens
            idx = nidx(bidx);               %% node indices for gens
        end
    end     %% methods
end         %% classdef
