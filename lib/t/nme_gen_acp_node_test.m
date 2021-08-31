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
        function nidxs = node_indices(obj, nm, dm, varargin)
            dme = obj.data_model_element(dm);
            b = dme.bus(dme.on);        %% gen bus index vector
            bt = dme.bus_etv(dme.on);   %% gen bus element type vector
            bidx = zeros(size(b));
            idx = zeros(size(b));
            for k = 1:dme.nbet
                if dm.elements.is_index_name(dme.bus_elm_types{k})
                    bus_dme = dm.elements.(dme.bus_elm_types{k});
                    beti = bus_dme.bus_eti;
                    nidx = nm.get_node_idx(dme.bus_elm_types{k});   %% node indices for bus element type
                    bidx(bt == beti) = bus_dme.i2on(b(bt == beti)); %% online bus indices for gens
                    idx(bt == beti) = nidx(bidx(bt == beti));   %% node indices for gens
                end
            end
            nidxs = {idx};
        end
    end     %% methods
end         %% classdef
