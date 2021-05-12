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
            f = dme.fbus(dme.on);
            t = dme.tbus(dme.on);
            bf = dme.fbus_type(dme.on);
            bt = dme.tbus_type(dme.on);
            fbidx = zeros(size(f));
            tbidx = zeros(size(t));
            fidx = zeros(size(f));
            tidx = zeros(size(t));

            bus_nld_dme = dm.elm_by_name('bus_nld');
            bus_ld_dme  = dm.elm_by_name('bus_ld');
            if ~isempty(bus_nld_dme)
                nidx_nld = nm.get_node_idx('bus_nld');  %% node indices for 'bus_nld'
                fbidx(bf == 1) = bus_nld_dme.i2on(f(bf == 1));  %% online bus_nld indices branch "from"
                tbidx(bt == 1) = bus_nld_dme.i2on(t(bt == 1));  %% online bus_nld indices branch "to"
                fidx(bf == 1) = nidx_nld(fbidx(bf == 1));   %% branch "from" port node indices
                tidx(bt == 1) = nidx_nld(tbidx(bt == 1));   %% branch "to" port node indices
            end
            if ~isempty(bus_ld_dme)
                nidx_ld  = nm.get_node_idx('bus_ld');   %% node indices for 'bus_ld'
                fbidx(bf == 2) = bus_ld_dme.i2on( f(bf == 2));  %% online bus_ld indices branch "from"
                tbidx(bt == 2) = bus_ld_dme.i2on( t(bt == 2));  %% online bus_ld indices branch "to"
                fidx(bf == 2) = nidx_ld( fbidx(bf == 2));   %% branch "from" port node indices
                tidx(bt == 2) = nidx_ld( tbidx(bt == 2));   %% branch "to" port node indices
            end
        end
    end     %% methods
end         %% classdef
