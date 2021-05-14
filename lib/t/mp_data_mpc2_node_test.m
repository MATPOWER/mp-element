classdef mp_data_mpc2_node_test < mp_data_mpc2
%MP_DATA_MPC2_NODE_TEST  Implementation of MATPOWER data model for MATPOWER case format v2

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
        function obj = mp_data_mpc2_node_test()
            %% call parent constructor
            obj@mp_data_mpc2();
            obj.element_classes = ...
                { @dme_bus_nld_mpc2_node_test, @dme_bus_ld_mpc2_node_test, @dme_gen_mpc2_node_test, ...
                    @dme_branch_mpc2_node_test };
        end

        %%-----  OPF methods  -----
        function [A, l, u, i] = branch_angle_diff_constraint(obj, ignore);
            branch_dme = obj.elm_by_name('branch');
            branch = branch_dme.get_table(obj);
            nb = 0;
            for k = 1:branch_dme.nbet
                bus_dme = obj.elm_by_name(branch_dme.bus_elm_types{k});
                if ~isempty(bus_dme)
                    nb = nb + bus_dme.n;
                end
            end
            mpopt = struct('opf', struct('ignore_angle_lim', ignore));
            [A, l, u, i]  = makeAang(obj.baseMVA, branch, nb, mpopt);
            if length(i)
                warning('OPF branch angle difference limits not implemented for this case.');
            end
        end
    end     %% methods
end         %% classdef
