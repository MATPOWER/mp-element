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
                bc = bus_nld_dme.bus_class;
                nidx_nld = nm.get_node_idx('bus_nld');  %% node indices for 'bus_nld'
                fbidx(bf == bc) = bus_nld_dme.i2on(f(bf == bc));%% online bus_nld indices branch "from"
                tbidx(bt == bc) = bus_nld_dme.i2on(t(bt == bc));%% online bus_nld indices branch "to"
                fidx(bf == bc) = nidx_nld(fbidx(bf == bc)); %% branch "from" port node indices
                tidx(bt == bc) = nidx_nld(tbidx(bt == bc)); %% branch "to" port node indices
            end
            if ~isempty(bus_ld_dme)
                bc = bus_ld_dme.bus_class;
                nidx_ld  = nm.get_node_idx('bus_ld');   %% node indices for 'bus_ld'
                fbidx(bf == bc) = bus_ld_dme.i2on( f(bf == bc));%% online bus_ld indices branch "from"
                tbidx(bt == bc) = bus_ld_dme.i2on( t(bt == bc));%% online bus_ld indices branch "to"
                fidx(bf == bc) = nidx_ld( fbidx(bf == bc)); %% branch "from" port node indices
                tidx(bt == bc) = nidx_ld( tbidx(bt == bc)); %% branch "to" port node indices
            end
        end
    end     %% methods
end         %% classdef
