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
            obj.element_classes = { @mpe_bus_dc, @mpe_gen_dc, @mpe_load_dc, @mpe_branch_dc, @mpe_shunt_dc };
            if isempty(obj.node)    %% skip if constructed from existing object
                obj.init_set_types();   %% should be called in mp_idx_manager
                                        %% constructor, if not for:
            end                         %% https://savannah.gnu.org/bugs/?52614
        end

        function obj = def_set_types(obj)
            def_set_types@mpe_network(obj);     %% call parent first
            obj.set_types.va = 'VOLTAGE VARS (va)';
            obj.set_types.z  = 'NON-VOLTAGE VARS (z)';
        end

        function obj = build_params(obj, nm, mpc)
            %% call parent to build individual element parameters
            build_params@mpe_network(obj, nm, mpc);

            %% aggregate parameters from individual elements
            obj.B = obj.stack_matrix_params('B', 1);
            obj.K = obj.stack_matrix_params('K', 0);
            obj.p = obj.stack_vector_params('p');
        end


        %%-----  PF methods  -----
        function ad = power_flow_aux_data(obj, mpc, mpopt)
            %% get model variables
            vvars = obj.model_vvars();
            zvars = obj.model_zvars();
            va = obj.params_var(vvars{1});
            z = obj.params_var(zvars{1});

            %% define constants
            [PQ, PV, REF, NONE] = idx_bus;

            %% get node types
            ntv = obj.power_flow_node_types(obj, mpc);
            ref = find(ntv == REF);     %% reference node indices
            pv  = find(ntv == PV );     %% PV node indices
            pq  = find(ntv == PQ );     %% PQ node indices

            %% get parameters
            [B, K, p] = obj.get_params();
            BB = obj.C * B * obj.C';
            pbus = -(obj.C * K * obj.D' * z + obj.C * p);
            branch_mpe = obj.mpe_by_name('branch');
            [Bf, pf] = branch_mpe.get_params(1:branch_mpe.nk, {'B', 'p'});

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
                'Bf',   Bf * branch_mpe.C', ...
                'Pbus', pbus, ...
                'Pfinj',pf  ...
            );
        end

        function [va, success, i, ad] = solve_power_flow(obj, mpc, mpopt)
            %% constant
            va_threshold = 1e5;     %% arbitrary threshold on |va| for declaring failure

            %% set up to trap non-singular matrix warnings
            [lastmsg, lastid] = lastwarn;
            lastwarn('');

            %% call parent
            [va, success, i, ad] = solve_power_flow@mpe_network(obj, mpc, mpopt);

            [msg, id] = lastwarn;
            %% Octave is not consistent in assigning proper warning id, so we'll just
            %% check for presence of *any* warning
            if ~isempty(msg) || max(abs(va)) > va_threshold
                success = 0;
            end

            %% restore warning state
            lastwarn(lastmsg, lastid);

            i = 1;                  %% not iterative
        end

        function add_pf_vars(obj, nm, om, mpc, mpopt)
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

        function add_pf_node_balance_constraints(obj, om, mpc, mpopt)
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
            branch_mpe = obj.mpe_by_name('branch');
            [Bbr, pbr] = branch_mpe.get_params(1:branch_mpe.nk, {'B', 'p'});
            om.userdata.Bf = Bbr * branch_mpe.C';
            om.userdata.Pfinj = pbr;
        end

        function names = opf_legacy_user_var_names(obj)
            names = {'Va', 'Pg'};
        end
    end     %% methods
end         %% classdef
