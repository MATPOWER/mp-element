classdef mme_bus_ld_pf_acps_node_test < mme_bus_nld_pf_acps_node_test

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        %% constructor
        function obj = mme_bus_ld_pf_acps_node_test()
            obj@mme_bus_nld_pf_acps_node_test();
            obj.name = 'bus_ld';
        end
    end     %% methods
end         %% classdef