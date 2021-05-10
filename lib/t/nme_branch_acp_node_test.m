classdef nme_branch_acp_node_test < nme_branch_acp

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function [fidx, tidx] = node_indices(obj, nm, dm, dme)
            bus_dme = dm.elm_by_name('bus_ld');
            nidx = nm.get_node_idx('bus_ld');   %% node indices for 'bus'
            f = dme.fbus(dme.on);
            t = dme.tbus(dme.on);
            fbidx = bus_dme.i2on(f);        %% online bus indices branch "from"
            tbidx = bus_dme.i2on(t);        %% online bus indices branch "to"
            fidx = nidx(fbidx);             %% branch "from" port node indices
            tidx = nidx(tbidx);             %% branch "to" port node indices
        end
    end     %% methods
end         %% classdef
