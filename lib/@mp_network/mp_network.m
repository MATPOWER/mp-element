classdef mpe_network < nm_element & mpe_container & mp_idx_manager% & mp_model
%MPE_NETWORK  Abstract base class for MATPOWER network model
%   Explicitly a subclass of NM_ELEMENT, MP_IDX_MANAGER and MPE_CONTAINER,
%   and implicitly assumed to be a subclass of MP_MODEL as well.

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        nme_port_map = [];      %% nme_port_map(k, 1:2), indices of 1st & last port for element k
        nme_z_map = [];         %% nme_z_map(k, 1:2), indices of 1st & last z var for element k
        nv = 0;                 %% total number of (real) v variables
        node = [];
        state = [];
    end
    
    methods
        %% constructor
        function obj = mpe_network()
            obj@nm_element();
            obj.name = 'network';
            obj.np = 0;     %% unknown number of ports at this point, init to 0
            obj.nk = 1;
            obj.nz = 0;     %% unknown number of z_ vars at this point, init to 0

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

        function obj = create_model(obj, dm, mpopt)
            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  Initialize MP_IDX_MANAGER fields here, if needed,
            %%              after object construction, but before object use.
            if isempty(obj.node)        %% only if not already initialized
                obj.init_set_types();
            end

            %%-----  HACK ALERT  -----
            %% This is a hack to deal with experimental
            %% mpopt.exp.sys_wide_zip_loads.pw/qw. MPOPT should be removed
            %% completely as an argument to create_model() once data and
            %% options are properly separated.
            if nargin == 3 && isfield(mpopt, 'exp') && ...
                    ~isempty(mpopt.exp) && ...
                    isfield(mpopt.exp, 'sys_wide_zip_loads') && ...
                    (~isempty(mpopt.exp.sys_wide_zip_loads.pw) || ...
                     ~isempty(mpopt.exp.sys_wide_zip_loads.qw))
                dm.mpc.sys_wide_zip_loads = mpopt.exp.sys_wide_zip_loads;
            end
            %%-----  end of HACK  -----

            %% create element objects for each class with data
            i = 0;
            for c = obj.element_classes
                nme = c{1}();       %% element constructor
                if nme.count(dm)
                    i = i + 1;
                    obj.elm_list{i} = nme;
                    obj.elm_map.(nme.name) = i;
                    obj.np = obj.np + nme.np * nme.nk;  %% number of ports
                    obj.nz = obj.nz + nme.nz * nme.nk;  %% number of z_ vars
                end
            end
            
            if obj.np ~= 0      %% skip for empty model
                %% create nodes and node voltage state variables
                obj.add_nodes(obj, dm);
            
                %% create non-voltage states and corresponding state variables
                obj.add_states(obj, dm);
            
                %% build params
                obj.build_params(obj, dm);
            end
        end

        function obj = add_nodes(obj, nm, dm)
            %% each element adds its nodes
            for nme = obj.elm_list
                nme{1}.add_nodes(obj, dm);
            end
            
            %% add voltage variables for each node
            obj.add_vvars(obj, dm);
        end

        function obj = add_states(obj, nm, dm)
            %% each element adds its states
            for nme = obj.elm_list
                nme{1}.add_states(obj, dm);
            end
            
            %% add state variables for each node
            obj.add_zvars(obj, dm);
        end

        function obj = build_params(obj, nm, dm)
            %% each element builds parameters, aggregate incidence matrices
            C = {};
            D = {};
            %% initialize nme_port_map, nme_z_map
            obj.nme_port_map = zeros(length(obj.elm_list), 2);
            obj.nme_z_map    = zeros(length(obj.elm_list), 2);
            pk = 1;     %% port counter
            zk = 1;     %% z-var counter
            for k = 1:length(obj.elm_list)
                nme = obj.elm_list{k};
                obj.nme_port_map(k, 1) = pk;        %% starting port index
                obj.nme_z_map(k, 1) = zk;           %% starting z-var index
                nme.build_params(obj, dm);
                C = horzcat(C, {nme.C});
                D = horzcat(D, {nme.D});
                pk = pk + nme.np * nme.nk;          %% increment port counter
                zk = zk + nme.nz * nme.nk;          %% increment z-var counter
                obj.nme_port_map(k, 2) = pk - 1;    %% ending port index
                obj.nme_z_map(k, 2)    = zk - 1;    %% ending z-var index
            end
            obj.C = horzcat(C{:});
            obj.D = horzcat(D{:});
        end

        function M = stack_matrix_params(obj, name, vnotz)
            nn = obj.getN('node');
            if vnotz
                nc = nn;
            else
                nc = obj.nz;
            end
            ii = {};
            jj = {};
            ss = {};
            last_i = 0;
            last_j = 0;
            for nme = obj.elm_list
                Mk = nme{1}.(name);
                if ~isempty(Mk)
                    [i, j, s] = find(Mk);
                    ii = horzcat(ii, i + last_i);
                    jj = horzcat(jj, j + last_j);
                    ss = horzcat(ss, s);
                end
                m = nme{1}.nk * nme{1}.np;      %% total number of ports for class
                if vnotz
                    n = m;
                else
                    n = nme{1}.nk * nme{1}.nz;  %% total number of states for class
                end
                last_i = last_i + m;
                last_j = last_j + n;
            end
            
            M = sparse(vertcat(ii{:}), vertcat(jj{:}), vertcat(ss{:}), last_i, last_j);
        end

        function v = stack_vector_params(obj, name)
            nn = obj.getN('node');
            vv = {};

            for nme = obj.elm_list
                vk = nme{1}.(name);
                if isempty(vk)
                    vk = zeros(nme{1}.nk * nme{1}.np, 1);
                end
                vv = horzcat(vv, vk);
            end
            v = vertcat(vv{:});
        end

        function obj = add_vvars(obj, nm, dm, idx)
            for k = 1:length(obj.node.order)
                nme = obj.elm_by_name(obj.node.order(k).name);
                nme.add_vvars(obj, dm, obj.state.order(k).idx);
            end
            for vtype = obj.model_vvars
                obj.nv = obj.nv + obj.getN(vtype{1});
            end
        end

        function obj = add_zvars(obj, nm, dm, idx)
            for k = 1:length(obj.state.order)
                nme = obj.elm_by_name(obj.state.order(k).name);
                nme.add_zvars(obj, dm, obj.state.order(k).idx);
            end
        end

        %%-----  mp_idx_manager methods  -----
        display(obj)
        
        [v0, vl, vu, vt] = params_var(obj, vtype, name, idx)
        
        function obj = def_set_types(obj)
            obj.set_types = struct(...
                    'node', 'NODES', ...
                    'state', 'STATES' ...
                );
        end

        function obj = init_set_types(obj)
            %% call parent to create base data structures for each type
            init_set_types@mp_idx_manager(obj);
            
            %% finish initializing data structures for each type
            es = struct();                          %% empty struct
                                                    %% variable types
            for vtype = horzcat(obj.model_vvars(), obj.model_zvars())
                assert(isfield(obj.set_types, vtype{1}), ...
                    'var type ''%'' is missing from def_set_types()', vtype{1});
                obj.(vtype{1}).data = struct( ...
                    'v0', es, ...
                    'vl', es, ...
                    'vu', es, ...
                    'vt', es );
            end
        end

        function add_node(obj, name, idx, N)
            %   obj.add_node(name, N)
            %   obj.add_node(name, idx_list, N)
            if ~iscell(idx)
                N = idx;
                idx = {};
            end

            %% add the named node set
            obj.add_named_set('node', name, idx, N);
        end

        function add_state(obj, name, idx, N)
            %   obj.add_state(name, N)
            %   obj.add_state(name, idx_list, N)
            if ~iscell(idx)
                N = idx;
                idx = {};
            end

            %% add the named state set
            obj.add_named_set('state', name, idx, N);
        end

        function add_var(obj, vtype, name, idx, varargin)
            %   obj.add_var(vtype, name, N, v0, vl, vu, vt)
            %   obj.add_var(vtype, name, N, v0, vl, vu)
            %   obj.add_var(vtype, name, N, v0, vl)
            %   obj.add_var(vtype, name, N, v0)
            %   obj.add_var(vtype, name, N)
            %   obj.add_var(vtype, name, idx_list, N, v0, vl, vu, vt)
            %   obj.add_var(vtype, name, idx_list, N, v0, vl, vu)
            %   obj.add_var(vtype, name, idx_list, N, v0, vl)
            %   obj.add_var(vtype, name, idx_list, N, v0)
            %   obj.add_var(vtype, name, idx_list, N)
            if iscell(idx)
                N = varargin{1};
                args = varargin(2:end);
            else
                N = idx;
                idx = {};
                args = varargin;
            end
            nargs = length(args);
    
            v0 = []; vl = []; vu = []; vt = [];
            if nargs >= 1
                v0 = args{1};
                if N > 1 && length(v0) == 1         %% expand from scalar as needed
                    v0 = v0 * ones(N, 1);
                end
                if nargs >= 2
                    vl = args{2};
                    if N > 1 && length(vl) == 1     %% expand from scalar as needed
                        vl = vl * ones(N, 1);
                    end
                    if nargs >= 3
                        vu = args{3};
                        if N > 1 && length(vu) == 1 %% expand from scalar as needed
                            vu = vu * ones(N, 1);
                        end
                        if nargs >= 4
                            vt = args{4};
                        end
                    end
                end
            end
            if isempty(v0)
                v0 = zeros(N, 1);   %% init to zero by default
            end
            if isempty(vl)
                vl = -Inf(N, 1);    %% unbounded below by default
            end
            if isempty(vu)
                vu = Inf(N, 1);     %% unbounded above by default
            end
            if isempty(vt) && N > 0
                vt = 'C';           %% all continuous by default
            end

            %% add the named variable set
            obj.add_named_set(vtype, name, idx, N);
            
            %% add type-specific data for var (v0, vl, vu, vt)
            d = obj.(vtype).data;
            obj.(vtype).data = [];
            if isempty(idx)
                d.v0.(name) = v0;       %% initial value
                d.vl.(name) = vl;       %% lower bound
                d.vu.(name) = vu;       %% upper bound
                d.vt.(name) = vt;       %% variable type
            else
                %% calls to substruct() are relatively expensive, so we
                %% pre-build the struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                d.v0 = subsasgn(d.v0, sc, v0);      %% initial value
                d.vl = subsasgn(d.vl, sc, vl);      %% lower bound
                d.vu = subsasgn(d.vu, sc, vu);      %% upper bound
                d.vt = subsasgn(d.vt, sc, vt);      %% variable type
            end
            obj.(vtype).data = d;
        end

        function varargout = get_node_idx(obj, name)
            %% [i1 iN] = obj.get_node_idx(name)
            %% nidx = obj.get_node_idx(name), where
            %%      nidx = [i1:iN]' or
            %%      nidx = {[i1(1):iN(1)]', ..., [i1(n):iN(n)]'}
            [varargout{1:nargout}] = obj.get_node_state_idx('node', name);
        end

        function varargout = get_state_idx(obj, name)
            %% [i1 iN] = obj.get_state_idx(name)
            %% nidx = obj.get_state_idx(name), where
            %%      nidx = [i1:iN]' or
            %%      nidx = {[i1(1):iN(1)]', ..., [i1(n):iN(n)]'}
            [varargout{1:nargout}] = obj.get_node_state_idx('state', name);
        end

        function [i1, iN] = get_node_state_idx(obj, node_state, name)
            %% private method
            %% [i1 iN] = obj.get_node_state_idx(node_state, name)
            %% nidx = obj.get_node_state_idx(node_state, name), where
            %%      nidx = [i1:iN]' or
            %%      nidx = {[i1(1):iN(1)]', ..., [i1(n):iN(n)]'}
            idx = obj.get_idx(node_state);
            i1 = idx.i1.(name);
            iN = idx.iN.(name);
            if nargout == 1
                N = length(i1);
                if N == 1
                    i1 = (i1:iN)';
                else
                    t = cell(1, N);
                    for k = 1:N
                        t{k} = (i1(k):iN(k))';
                    end
                    i1 = t;
                end
            end
        end

        function [ref, pv, pq] = node_types(obj, nm, dm)
            %%           ntv = node_types(obj, nm, dm)
            %% [ref, pv, pq] = node_types(obj, nm, dm)
            %% get node type vector from each node-creating NME
            tt = cell(length(obj.node.order), 1);
            for k = 1:length(obj.node.order)
                nme = obj.elm_by_name(obj.node.order(k).name);
                tt{k} = nme.node_types(obj, dm);
            end
            ntv = vertcat(tt{:});       %% concatenate into a single vector

            if nargout > 1
                ref = dm.node_type_ref(ntv);    %% reference node indices
                pv  = dm.node_type_pv(ntv);     %% PV node indices
                pq  = dm.node_type_pq(ntv);     %% PQ node indices
            else
                ref = ntv;
            end
        end

        %%-----  PF methods  -----
        function add_pf_constraints(obj, nm, om, dm, mpopt)
            %% system constraints
            obj.add_pf_system_constraints(om, dm, mpopt);
            
%             %% each element adds its PF constraints
%             for nme = obj.elm_list
%                 nme{1}.add_pf_constraints(nm, om, dm, mpopt);
%             end
        end

        function add_pf_system_constraints(obj, om, dm, mpopt)
            %% can be overridden to add additional system constraints

            %% node balance constraints
            obj.add_pf_node_balance_constraints(om, dm, mpopt);
        end

        opt = solve_opts_power_flow(obj, om, dm, mpopt)


        %%-----  OPF methods  -----
        om = setup_opf(obj, dm, mpopt)
        
        [x, success, i] = solve_opf(obj, dm, mpopt)
        
        function add_opf_vars(obj, nm, om, dm, mpopt)
            vars = horzcat(obj.model_vvars(), obj.model_zvars());
            for vtype = vars
                st = obj.(vtype{1});    %% set type
                for k = 1:st.NS
                    name = st.order(k).name;
                    if isempty(st.order(k).idx)
                        d = st.data;
                        om.add_var(name, st.idx.N.(name), d.v0.(name), d.vl.(name), d.vu.(name), d.vt.(name));
                    else
                        error('handling of indexed sets not implmented here (yet)');
                    end
                end
            end
            
            %% each element adds its OPF variables
            for nme = obj.elm_list
                nme{1}.add_opf_vars(nm, om, dm, mpopt);
            end
            
            %% legacy user-defined variables
            obj.add_opf_legacy_user_vars(om, dm, mpopt);
        end

        function add_opf_constraints(obj, nm, om, dm, mpopt)
            %% system constraints
            obj.add_opf_system_constraints(om, dm, mpopt);
            
            %% each element adds its OPF constraints
            for nme = obj.elm_list
                nme{1}.add_opf_constraints(nm, om, dm, mpopt);
            end
        end

        function add_opf_costs(obj, nm, om, dm, mpopt)
            %% system costs
            obj.add_opf_system_costs(om, dm, mpopt);
            
            %% each element adds its OPF costs
            for nme = obj.elm_list
                nme{1}.add_opf_costs(nm, om, dm, mpopt);
            end
        end

        function add_opf_system_constraints(obj, om, dm, mpopt)
            %% can be overridden to add additional system constraints

            %% node balance constraints
            obj.add_opf_node_balance_constraints(om);

            %% legacy user-defined constraints
            obj.add_opf_legacy_user_constraints(om, dm, mpopt);
        end

        function add_opf_system_costs(obj, om, dm, mpopt)
            %% can be overridden to add additional system costs

            %% legacy user-defined costs
            obj.add_opf_legacy_user_costs(om, mpopt);
        end

        function add_opf_legacy_user_vars(obj, om, dm, mpopt)
            uv_names = obj.opf_legacy_user_var_names();
            nx = sum(cellfun(@(x)om.getN('var', x), uv_names));
            [uv, z] = dm.opf_legacy_user_vars(uv_names, nx, mpopt);

            %% save data
            om.userdata = nested_struct_copy(om.userdata, uv);

            %% add any user-defined vars
            if z.nz > 0
                om.add_var('z', z.nz, z.z0, z.zl, z.zu);
                om.userdata.user_vars{end+1} = 'z';
            end
        end

        function add_opf_legacy_user_constraints(obj, om, dm, mpopt)
            %% user-defined linear constraints
            if om.userdata.nlin
                om.add_lin_constraint('usr', om.userdata.A, om.userdata.l, ...
                    om.userdata.u, om.userdata.user_vars);
                om.userdata = rmfield(om.userdata, {'A', 'l', 'u', 'nlin'});
            end
        end

        function add_opf_legacy_user_costs(obj, om, mpopt)
            if om.userdata.nw
                user_cost.N = om.userdata.N;
                user_cost.Cw = om.userdata.Cw;
                if ~isempty(om.userdata.fparm)
                    user_cost.dd = om.userdata.fparm(:, 1);
                    user_cost.rh = om.userdata.fparm(:, 2);
                    user_cost.kk = om.userdata.fparm(:, 3);
                    user_cost.mm = om.userdata.fparm(:, 4);
                end
                if ~isempty(om.userdata.H)
                    user_cost.H = om.userdata.H;
                end
                om.add_legacy_cost('usr', user_cost, om.userdata.user_vars);
                om.userdata = rmfield(om.userdata, {'N', 'fparm', 'H', 'Cw', 'nw'});
            end

            %% implement legacy user costs using quadratic or general non-linear costs
            dc = strcmp(upper(mpopt.model), 'DC');
            cp = om.params_legacy_cost();   %% construct/fetch the parameters
            [N, H, Cw, rh, mm] = deal(cp.N, cp.H, cp.Cw, cp.rh, cp.mm);
            [nw, nx] = size(N);
            if nw
                if any(cp.dd ~= 1) || any(cp.kk)    %% not simple quadratic form
                    if dc                           %% (includes "dead zone" or
                        if any(cp.dd ~= 1)          %%  quadratic "penalty")
                            error('mpe_network/add_opf_legacy_user_costs: DC OPF can only handle legacy user-defined costs with d = 1');
                        end
                        if any(cp.kk)
                            error('mpe_network/add_opf_legacy_user_costs: DC OPF can only handle legacy user-defined costs with no "dead zone", i.e. k = 0');
                        end
                    else
                        %% use general nonlinear cost to implement legacy user cost
                        legacy_cost_fcn = @(x)opf_legacy_user_cost_fcn(x, cp);
                        om.add_nln_cost('usr', 1, legacy_cost_fcn);
                    end
                else                                %% simple quadratic form
                    %% use a quadratic cost to implement legacy user cost
                    %% f = 1/2 * w'*H*w + Cw'*w, where w = diag(mm)*(N*x - rh)
                    %% Let: MN = diag(mm)*N
                    %%      MR = M * rh
                    %%      HMR  = H  * MR;
                    %%      HtMR = H' * MR;
                    %%  =>   w = MN*x - MR
                    %% f = 1/2 * (MN*x - MR)'*H*(MN*x - MR) + Cw'*(MN*x - MR)
                    %%   = 1/2 * x'*MN'*H*MN*x +
                    %%          (Cw'*MN - 1/2 * MR'*(H+H')*MN)*x +
                    %%          1/2 * MR'*H*MR - Cw'*MR
                    %%   = 1/2 * x'*Q*w + c'*x + k

                    [N, H, Cw, rh, mm] = deal(cp.N, cp.H, cp.Cw, cp.rh, cp.mm);
                    nw = size(N, 1);            %% number of general cost vars, w
                    M    = sparse(1:nw, 1:nw, mm, nw, nw);
                    MN   = M * N;
                    MR   = M * rh;
                    HMR  = H  * MR;
                    HtMR = H' * MR;
                    Q = MN' * H * MN;
                    c = full(MN' * (Cw - 1/2*(HMR+HtMR)));
                    k = (1/2 * HtMR - Cw)' * MR;
                    om.add_quad_cost('usr', Q, c, k);
                end
            end
        end
    end     %% methods
end         %% classdef