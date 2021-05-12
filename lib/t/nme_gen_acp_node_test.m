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
            b = dme.bus(dme.on);        %% gen "bus"
            bt = dme.bus_type(dme.on);  %% gen "bus type"
            bidx = zeros(size(b));
            idx = zeros(size(b));
            bus_nld_dme = dm.elm_by_name('bus_nld');
            bus_ld_dme  = dm.elm_by_name('bus_ld');
            if ~isempty(bus_nld_dme)
                nidx_nld = nm.get_node_idx('bus_nld');  %% node indices for 'bus_nld'
                bidx(bt == 1) = bus_nld_dme.i2on(b(bt == 1));   %% online bus_nld indices for gens
                idx(bt == 1) = nidx_nld(bidx(bt == 1)); %% node indices for gens
            end
            if ~isempty(bus_ld_dme)
                nidx_ld  = nm.get_node_idx('bus_ld');   %% node indices for 'bus_ld'
                bidx(bt == 2) = bus_ld_dme.i2on( b(bt == 2));   %% online bus_ld indices for gens
                idx(bt == 2) = nidx_ld( bidx(bt == 2)); %% node indices for gens
            end
        end
    end     %% methods
end         %% classdef
