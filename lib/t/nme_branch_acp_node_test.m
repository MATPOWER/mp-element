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
            nidxs = node_indices@nme_branch_acp(obj, nm, dm, dme.cxn_type, dme.cxn_idx_prop, dme.cxn_type_prop);
        end
    end     %% methods
end         %% classdef
