classdef nme_branch_acp_node_test < nme_branch_acp

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        function nidxs = node_indices(obj, nm, dm, varargin)
            dme = obj.data_model_element(dm);
            f = dme.fbus(dme.on);
            t = dme.tbus(dme.on);
            bf = dme.fbus_etv(dme.on);
            bt = dme.tbus_etv(dme.on);
            fbidx = zeros(size(f));
            tbidx = zeros(size(t));
            fidx = zeros(size(f));
            tidx = zeros(size(t));

            for k = 1:dme.nbet
                if dm.elements.is_index_name(dme.bus_elm_types{k})
                    bus_dme = dm.elements.(dme.bus_elm_types{k});
                    beti = bus_dme.bus_eti;
                    nidx = nm.get_node_idx(dme.bus_elm_types{k});   %% node indices for bus element type
                    fbidx(bf == beti) = bus_dme.i2on(f(bf == beti));%% online bus indices of branch "from"
                    tbidx(bt == beti) = bus_dme.i2on(t(bt == beti));%% online bus indices of branch "to"
                    fidx(bf == beti) = nidx(fbidx(bf == beti)); %% branch "from" port node indices
                    tidx(bt == beti) = nidx(tbidx(bt == beti)); %% branch "to" port node indices
                end
            end
            nidxs = {fidx, tidx};
        end
    end     %% methods
end         %% classdef
