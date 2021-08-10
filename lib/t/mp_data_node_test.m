classdef mp_data_node_test < mp_data
%MP_DATA_NODE_TEST  Implementation of MATPOWER data model for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020-2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        %% constructor
        function obj = mp_data_node_test()
            %% call parent constructor
            obj@mp_data();
            obj.element_classes = ...
                { @dme_bus_nld_node_test, @dme_bus_ld_node_test, ...
                    @dme_gen_node_test, @dme_branch_node_test };
        end
    end     %% methods
end         %% classdef
