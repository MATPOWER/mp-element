classdef dc_aggregate < mp_aggregate & dc_model

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        va = [];
        z  = [];
    end
    
    methods
        %% constructor
        function obj = dc_aggregate()
            obj@mp_aggregate();
            obj.element_classes = { @dc_bus, @dc_gen, @dc_load, @dc_branch, @dc_shunt };
            if isempty(obj.node)    %% skip if constructed from existing object
                obj.init_set_types();   %% should be called in mp_idx_manager
                                        %% constructor, if not for:
            end                         %% https://savannah.gnu.org/bugs/?52614
        end

        function obj = def_set_types(obj)
            def_set_types@mp_aggregate(obj);        %% call parent first
            obj.set_types.va = 'VOLTAGE VARS (va)';
            obj.set_types.z  = 'NON-VOLTAGE VARS (z)';
        end

        function obj = build_params(obj, asm, mpc)
            %% call parent to build individual element parameters
            build_params@mp_aggregate(obj, asm, mpc);

            %% aggregate parameters from individual elements
            obj.B = obj.stack_matrix_params('B', 1);
            obj.K = obj.stack_matrix_params('K', 0);
            obj.p = obj.stack_vector_params('p');
        end

        function add_opf_node_balance_constraints(obj, om)
            [B, K, p] = obj.get_params();

            %% power balance constraints
            C = obj.getC();
            Amis = C * [B*C' K*obj.getD('tr')];
            bmis = -C * p;
            om.add_lin_constraint('Pmis', Amis, bmis, bmis, ...
                                {obj.va.order(:).name obj.z.order(:).name});

            %% user data
            branch_mpe = obj.mpe_by_name('branch');
            [Bbr, pbr] = branch_mpe.get_params(1:branch_mpe.nk, {'B', 'p'});
            Cbrt = branch_mpe.getC('tr');
            om.userdata.Bf = Bbr * Cbrt;
            om.userdata.Pfinj = pbr;
        end
    end     %% methods
end         %% classdef
