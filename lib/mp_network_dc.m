classdef mp_network_dc < mp_network & mp_form_dc

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
        function obj = mp_network_dc()
            obj@mp_network();
            obj.element_classes = { ...
                @nme_bus_dc, @nme_gen_dc, @nme_load_dc, ...
                    @nme_branch_dc, @nme_shunt_dc };

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  INIT_SET_TYPES() is called explicitly as needed
            %%              (if obj.node is empty) in BUILD() and DISPLAY(),
            %%              after object construction, but before object use.
        end

        function obj = def_set_types(obj)
            def_set_types@mp_network(obj);      %% call parent first
            obj.set_types.va = 'VOLTAGE VARS (va)';
            obj.set_types.z  = 'NON-VOLTAGE VARS (z)';
        end

        function obj = build_params(obj, nm, dm)
            %% call parent to build individual element parameters
            build_params@mp_network(obj, nm, dm);

            %% aggregate parameters from individual elements
            obj.B = obj.stack_matrix_params('B', 1);
            obj.K = obj.stack_matrix_params('K', 0);
            obj.p = obj.stack_vector_params('p');
        end

        function obj = port_inj_soln(obj)
            %% compute port injections
            obj.soln.gp = obj.port_inj_power(obj.soln.x);
        end


        %%-----  PF methods  -----
        function ad = pf_aux_data(obj, dm, mpopt)
            %% get model variables
            vvars = obj.model_vvars();
            zvars = obj.model_zvars();
            va = obj.params_var(vvars{1});
            z = obj.params_var(zvars{1});

            %% get node types
            [ref, pv, pq] = obj.node_types(obj, dm);

            %% get parameters
            [B, K, p] = obj.get_params();

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
                'B',    obj.C * B * obj.C', ...
                'Pbus', -(obj.C * K * obj.D' * z + obj.C * p) ...
            );
        end

        function opt = pf_solve_opts(obj, mm, dm, mpopt)
            %% TO DO: move pf.alg to pf.ac.solver and add a
            %%        pf.dc.solver to set the 'leq_opt.solver' option here
            opt = struct( ...
                'verbose',  mpopt.verbose, ...
                'leq_opt',  struct('thresh', 1e5)   );
        end

        function obj = pf_add_vars(obj, mm, nm, dm, mpopt)
            %% get model variables
            vvars = obj.model_vvars();

            %% index vectors
            ad = mm.get_userdata('aux_data');
            pvq = [ad.pv; ad.pq];

            %% voltage angles
            st = obj.(vvars{1});
            for k = 1:st.NS
                name = st.order(k).name;
                if isempty(st.order(k).idx)
                    d = st.data;
                    mm.add_var(name, ad.npv+ad.npq, d.v0.(name)(pvq), d.vl.(name)(pvq), d.vu.(name)(pvq));
                else
                    error('mp_network_dc/pf_add_vars: handling of indexed sets not implmented here (yet)');
                end
            end
        end

        function [vx, z, x] = pf_convert_x(obj, mmx, ad, only_v)
            %% x = obj.pf_convert(mmx, ad)
            %% [v, z] = obj.pf_convert(mmx, ad)
            %% [v, z, x] = obj.pf_convert(mmx, ad)
            %% ... = obj.pf_convert(mmx, ad, only_v)

            %% update v_, z_ from mmx
            vx = ad.va;
            vx([ad.pv; ad.pq]) = mmx(1:ad.npv+ad.npq);      %% va
            z = ad.z;

            %% update z, if requested
            if nargin < 4 || ~only_v
                z = obj.pf_update_z(vx, z, ad);
            end

            %% prepare return values
            if nargout < 2
                vx = [vx; z];
            elseif nargout > 2
                x = [vx; z];
            end
        end

        function z = pf_update_z(obj, v, z, ad)
            %% update/allocate slack node active power injections
            
            %% coefficient matrix for power injection states at slack bus
            CC = obj.C(ad.ref, :) * obj.get_params([], 'K') * obj.D';
            jr = find(any(CC, 1));  %% indices of corresponding states

            %% power injections at slack nodes
            idx = find(any(obj.C(ad.ref, :), 1));  %% ports connected to slack nodes
            Pref = obj.C(ad.ref, idx) * obj.port_inj_power([v; z], 1, idx);

            %% allocate active power at slack nodes to 1st direct inj state
            %% find all z (except first one) with direct injection at each
            %% slack node
            [i, j] = find(CC);
            if size(i, 2) > 1, i = i'; j = j'; end
            ij = sortrows([i j]);       %% 1st state comes 1st for each node
            [~, k1] = unique(ij(:, 1), 'first');%% index of 1st entry for each node
            %% all included states that are not 1st at their node
            jn = unique(ij(~ismember(1:length(i), k1), 2));

            %% if we have extra states (more than 1) for any node(s)
            if ~isempty(jn)
                %% augment update equation CC * (z - zprev) = -Pref with
                %% additional rows to force these states to remain fixed
                I = speye(obj.nz);
                CC = [CC; I(jn, :)];
                Pref = [Pref; zeros(length(jn), 1)];
            end

            %% update z for active injections at slack nodes
            z(jr) = z(jr) - CC(:, jr) \ Pref;
        end

        function obj = pf_add_node_balance_constraints(obj, mm, dm, mpopt)
            ad = mm.get_userdata('aux_data');
            pvq = [ad.pv; ad.pq];

            %% power balance constraints
            A = ad.B(pvq, pvq);
            b = (ad.Pbus(pvq) - ad.B(pvq, ad.ref) * ad.va(ad.ref));
            mm.add_lin_constraint('Pmis', A, b, b);
        end


        %%-----  OPF methods  -----
        function [vx, z, x] = opf_convert_x(obj, mmx, ad)
            %% convert (real) math model x to network model x
            if nargout < 2
                vx = mmx(1:obj.nv+obj.nz);
            else
                vx = mmx(1:obj.nv);
                z  = mmx(obj.nv+1:obj.nv+obj.nz);
                if nargout > 2
                    x = [vx; z];
                end
            end
        end

        function opf_add_node_balance_constraints(obj, mm)
            [B, K, p] = obj.get_params();

            %% power balance constraints
            C = obj.C;
            Amis = C * [B*C' K*obj.D'];
            bmis = -C * p;
            mm.add_lin_constraint('Pmis', Amis, bmis, bmis, ...
                                {obj.va.order(:).name obj.z.order(:).name});
        end

        function opf_add_system_costs(obj, mm, dm, mpopt)
            %% can be overridden to add additional system costs

            %% legacy user-defined costs
            obj.opf_add_legacy_user_costs(mm, dm, 1);
        end

        function names = opf_legacy_user_var_names(obj)
            names = {'Va', 'Pg'};
        end

        function opt = opf_solve_opts(obj, mm, dm, mpopt)
            opt = mpopt2qpopt(mpopt, mm.problem_type());

            switch opt.alg
                case {'MIPS', 'IPOPT'}
                    if mpopt.opf.start < 2
                        %% initialize interior point
                        x0 = obj.opf_interior_x0(mm, dm);

                        %% set voltages
                        %% Va equal to angle of 1st ref bus
                        vv = mm.get_idx();
                        bus_dme = dm.elements.bus;
                        Varefs = bus_dme.Va0(find(bus_dme.isref));
                        x0(vv.i1.Va:vv.iN.Va) = Varefs(1);  %% angles set to first reference angle

                        opt.x0 = x0;
                    end
                case 'OSQP'
                    opt.x0 = [];        %% disable provided starting point
            end
        end
    end     %% methods
end         %% classdef
