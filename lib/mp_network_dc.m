classdef mpe_network_dc < mpe_network & mp_model_dc

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
        function obj = mpe_network_dc()
            obj@mpe_network();
            obj.element_classes = { ...
                @nme_bus_dc, @nme_gen_dc, @nme_load_dc, ...
                    @nme_branch_dc, @nme_shunt_dc };

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  INIT_SET_TYPES() is called explicitly as needed
            %%              (if obj.node is empty) in CREATE_MODEL() and
            %%              DISPLAY(), after object construction, but before
            %%              object use.
        end

        function obj = def_set_types(obj)
            def_set_types@mpe_network(obj);     %% call parent first
            obj.set_types.va = 'VOLTAGE VARS (va)';
            obj.set_types.z  = 'NON-VOLTAGE VARS (z)';
        end

        function obj = build_params(obj, nm, dm)
            %% call parent to build individual element parameters
            build_params@mpe_network(obj, nm, dm);

            %% aggregate parameters from individual elements
            obj.B = obj.stack_matrix_params('B', 1);
            obj.K = obj.stack_matrix_params('K', 0);
            obj.p = obj.stack_vector_params('p');
        end


        %%-----  PF methods  -----
        function ad = power_flow_aux_data(obj, dm, mpopt)
            %% get model variables
            vvars = obj.model_vvars();
            zvars = obj.model_zvars();
            va = obj.params_var(vvars{1});
            z = obj.params_var(zvars{1});

            %% get node types
            [ref, pv, pq] = obj.node_types(obj, dm);

            %% get parameters
            [B, K, p] = obj.get_params();
            BB = obj.C * B * obj.C';
            pbus = -(obj.C * K * obj.D' * z + obj.C * p);
            branch_nme = obj.elm_by_name('branch');
            [Bf, pf] = branch_nme.get_params(1:branch_nme.nk, {'B', 'p'});

            %% create aux_data struct
            ad = struct( ...
                'va', va, ...               %% initial value of va
                'z', z, ...                 %% initial value of z
                'ref',  ref, ...            %% REF node indices
                'nref', length(ref), ...    %% number of REF nodes
                'pv',   pv, ...             %% PV node indices
                'npv',  length(pv), ...     %% number of PV nodes
                'pq',   pq, ...             %% PQ node indices
                'npq',  length(pq), ...     %% number of PQ nodes
                'B',    BB, ...
                'Bf',   Bf * branch_nme.C', ...
                'Pbus', pbus, ...
                'Pfinj',pf  ...
            );
        end

        function opt = solve_opts_power_flow(obj, om, dm, mpopt)
            %% TO DO: move pf.alg to pf.ac.solver and add a
            %%        pf.dc.solver to set the 'leq_opt.solver' option here
            opt = struct( ...
                'verbose',  mpopt.verbose, ...
                'leq_opt',  struct('thresh', 1e5)   );
        end

        function add_pf_vars(obj, nm, om, dm, mpopt)
            %% get model variables
            vvars = obj.model_vvars();

            %% index vectors
            ad = om.get_userdata('power_flow_aux_data');
            pvq = [ad.pv; ad.pq];

            %% voltage angles
            st = obj.(vvars{1});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    om.add_var(name, ad.npv+ad.npq, d.v0.(name)(pvq), d.vl.(name)(pvq), d.vu.(name)(pvq));
                else
                    error('handling of indexed sets not implmented here (yet)');
                end
            end
        end

        function [va, z] = pfx2vz(obj, x, ad)
            %% update v_, z_ from x
            va = ad.va;
            va([ad.pv; ad.pq]) = x(1:ad.npv+ad.npq);        %% va
            z = ad.z;
        end

        function add_pf_node_balance_constraints(obj, om, dm, mpopt)
            ad = om.get_userdata('power_flow_aux_data');
            pvq = [ad.pv; ad.pq];

            %% power balance constraints
            A = ad.B(pvq, pvq);
            b = (ad.Pbus(pvq) - ad.B(pvq, ad.ref) * ad.va(ad.ref));
            om.add_lin_constraint('Pmis', A, b, b);
        end


        %%-----  OPF methods  -----
        function add_opf_node_balance_constraints(obj, om)
            [B, K, p] = obj.get_params();

            %% power balance constraints
            C = obj.C;
            Amis = C * [B*C' K*obj.D'];
            bmis = -C * p;
            om.add_lin_constraint('Pmis', Amis, bmis, bmis, ...
                                {obj.va.order(:).name obj.z.order(:).name});

            %% user data
            branch_nme = obj.elm_by_name('branch');
            [Bbr, pbr] = branch_nme.get_params(1:branch_nme.nk, {'B', 'p'});
            om.userdata.Bf = Bbr * branch_nme.C';
            om.userdata.Pfinj = pbr;
        end

        function names = opf_legacy_user_var_names(obj)
            names = {'Va', 'Pg'};
        end
    end     %% methods
end         %% classdef