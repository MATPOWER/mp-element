classdef mpe_network < mp_element & mp_idx_manager% & mp_model
%MPE_NETWORK  Abstract class, explicitly a subclass of MP_ELEMENT and
%             MP_IDX_MANAGER and implicitly assumed to be a subclass of
%             MP_MODEL as well

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        element_classes = {};   %% classes for individual element types
                                %% filled in by subclass constructor
        mpe_list = {};          %% cell array of mp_element objects
        mpe_map  = struct();    %% key = element name, val = index into mpe_list
        mpe_port_map = [];      %% mpe_port_map(k, 1:2), indices of 1st & last port for element k
        mpe_z_map = [];         %% mpe_z_map(k, 1:2), indices of 1st & last z var for element k
        nv = 0;                 %% total number of (real) v variables
        node = [];
        state = [];
    end
    
    methods
        %% constructor
        function obj = mpe_network()
            obj@mp_element();
            obj.name = 'network';
            obj.dm_table = '';
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

        function obj = modify_element_classes(obj, class_list)
            %% each element in class_list is either:
            %%  1 - a handle to a constructor to be appended to
            %%      obj.element_classes, or
            %%  2 - a 2-element cell array {A,B} where A is a handle to
            %%      a constructor to replace any element E in the list for
            %%      which isa(E(), B) is true, i.e. B is a char array
            if ~iscell(class_list)
                class_list = {class_list};
            end
            ec = obj.element_classes;   %% list to be updated
            ec0 = {ec{:}};              %% unmodified copy of original list
            for k = 1:length(class_list)
                c = class_list{k};
                if iscell(c)        %% it's a 2-d cell array
                    i = find(cellfun(@(e)isa(e(), c{2}), ec0)); %% find c{2}
                    if ~isempty(i)
                        ec{i} = c{1};                   %% replace with c{1}
                    end
                else                %% it's a single function handle
                    ec{end+1} = c;  %%      append it
                end
            end
            obj.element_classes = ec;
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
                mpe = c{1}();       %% element constructor
                if mpe.count(dm)
                    i = i + 1;
                    obj.mpe_list{i} = mpe;
                    obj.mpe_map.(mpe.name) = i;
                    obj.np = obj.np + mpe.np * mpe.nk;  %% number of ports
                    obj.nz = obj.nz + mpe.nz * mpe.nk;  %% number of z_ vars
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

        function mpe = mpe_by_name(obj, name)
            mpe = obj.mpe_list{obj.mpe_map.(name)};
        end

        function obj = add_nodes(obj, nm, dm)
            %% each element adds its nodes
            for mpe = obj.mpe_list
                mpe{1}.add_nodes(obj, dm);
            end
            
            %% add voltage variables for each node
            obj.add_vvars(obj, dm);
        end

        function obj = add_states(obj, nm, dm)
            %% each element adds its states
            for mpe = obj.mpe_list
                mpe{1}.add_states(obj, dm);
            end
            
            %% add state variables for each node
            obj.add_zvars(obj, dm);
        end

        function obj = build_params(obj, nm, dm)
            %% each element builds parameters, aggregate incidence matrices
            C = {};
            D = {};
            %% initialize mpe_port_map, mpe_z_map
            obj.mpe_port_map = zeros(length(obj.mpe_list), 2);
            obj.mpe_z_map    = zeros(length(obj.mpe_list), 2);
            pk = 1;     %% port counter
            zk = 1;     %% z-var counter
            for k = 1:length(obj.mpe_list)
                mpe = obj.mpe_list{k};
                obj.mpe_port_map(k, 1) = pk;        %% starting port index
                obj.mpe_z_map(k, 1) = zk;           %% starting z-var index
                mpe.build_params(obj, dm);
                C = horzcat(C, {mpe.C});
                D = horzcat(D, {mpe.D});
                pk = pk + mpe.np * mpe.nk;          %% increment port counter
                zk = zk + mpe.nz * mpe.nk;          %% increment z-var counter
                obj.mpe_port_map(k, 2) = pk - 1;    %% ending port index
                obj.mpe_z_map(k, 2)    = zk - 1;    %% ending z-var index
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
            for mpe = obj.mpe_list
                Mk = mpe{1}.(name);
                if ~isempty(Mk)
                    [i, j, s] = find(Mk);
                    ii = horzcat(ii, i + last_i);
                    jj = horzcat(jj, j + last_j);
                    ss = horzcat(ss, s);
                end
                m = mpe{1}.nk * mpe{1}.np;      %% total number of ports for class
                if vnotz
                    n = m;
                else
                    n = mpe{1}.nk * mpe{1}.nz;  %% total number of states for class
                end
                last_i = last_i + m;
                last_j = last_j + n;
            end
            
            M = sparse(vertcat(ii{:}), vertcat(jj{:}), vertcat(ss{:}), last_i, last_j);
        end

        function v = stack_vector_params(obj, name)
            nn = obj.getN('node');
            vv = {};

            for mpe = obj.mpe_list
                vk = mpe{1}.(name);
                if isempty(vk)
                    vk = zeros(mpe{1}.nk * mpe{1}.np, 1);
                end
                vv = horzcat(vv, vk);
            end
            v = vertcat(vv{:});
        end

        function obj = add_vvars(obj, nm, dm, idx)
            for k = 1:length(obj.node.order)
                mpe = obj.mpe_by_name(obj.node.order(k).name);
                mpe.add_vvars(obj, dm, obj.state.order(k).idx);
            end
            for vtype = obj.model_vvars
                obj.nv = obj.nv + obj.getN(vtype{1});
            end
        end

        function obj = add_zvars(obj, nm, dm, idx)
            for k = 1:length(obj.state.order)
                mpe = obj.mpe_by_name(obj.state.order(k).name);
                mpe.add_zvars(obj, dm, obj.state.order(k).idx);
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
            obj.node.data = struct( ...             %% node type
                'idx2ID', es, ...
                'ID2idx', es );

            obj.state.data = struct( ...            %% state type
                'idx2ID', es, ...
                'ID2idx', es );
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

        function add_node(obj, name, idx, varargin)
            %   obj.add_node(name, N, IDs)
            %   obj.add_node(name, N)
            %   obj.add_node(name, idx_list, N, IDs)
            %   obj.add_node(name, idx_list, N)
            if iscell(idx)
                N = varargin{1};
                args = varargin(2:end);
            else
                N = idx;
                idx = {};
                args = varargin;
            end
            nargs = length(args);
            idx2ID = [];
            if nargs >= 1;
                idx2ID = args{1};
                if length(idx2ID) ~= N
                    error('mpe_network/add_node: length of IDs vector (%d) must equal N (%d) ', length(idx2ID), N);
                end
            end
            if isempty(idx2ID)
                idx2ID = [1:N]';
            end

            %% create reverse mapping
            ID2idx = sparse(idx2ID, ones(N, 1), 1:N, max(idx2ID), 1);

            %% add the named node set
            obj.add_named_set('node', name, idx, N);
            
            %% add type-specific data for nodes (idx2ID, ID2idx)
            if isempty(idx)
                obj.node.data.idx2ID.(name) = idx2ID;
                obj.node.data.ID2idx.(name) = ID2idx;
            else
                %% calls to substruct() are relatively expensive, so we
                %% pre-build the struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                obj.node.data.idx2ID = subsasgn(obj.node.data.idx2ID, sc, idx2ID);
                obj.node.data.ID2idx = subsasgn(obj.node.data.ID2idx, sc, ID2idx);
            end
        end

        function add_state(obj, name, idx, varargin)
            %   obj.add_state(name, N, IDs)
            %   obj.add_state(name, N)
            %   obj.add_state(name, idx_list, N, IDs)
            %   obj.add_state(name, idx_list, N)
            if iscell(idx)
                N = varargin{1};
                args = varargin(2:end);
            else
                N = idx;
                idx = {};
                args = varargin;
            end
            nargs = length(args);
            idx2ID = [];
            if nargs >= 1;
                idx2ID = args{1};
                if length(idx2ID) ~= N
                    error('mpe_network/add_state: length of IDs vector (%d) must equal N (%d) ', length(idx2ID), N);
                end
            end
            if isempty(idx2ID)
                idx2ID = [1:N]';
            end

            %% create reverse mapping
            ID2idx = sparse(idx2ID, ones(N, 1), 1:N, max(idx2ID), 1);

            %% add the named state set
            obj.add_named_set('state', name, idx, N);
            
            %% add type-specific data for states (idx2ID, ID2idx)
            if isempty(idx)
                obj.state.data.idx2ID.(name) = idx2ID;
                obj.state.data.ID2idx.(name) = ID2idx;
            else
                %% calls to substruct() are relatively expensive, so we
                %% pre-build the struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                obj.state.data.idx2ID = subsasgn(obj.state.data.idx2ID, sc, idx2ID);
                obj.state.data.ID2idx = subsasgn(obj.state.data.ID2idx, sc, ID2idx);
            end
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


        %%-----  PF methods  -----
        function ntv = power_flow_node_types(obj, nm, dm, idx)
%         function [ntv, nts] = power_flow_node_types(obj, nm, dm, idx)
            %% create empty cell array for node type vectors
            tt = cell(length(obj.node.order), 1);
            
            %% get node type vector from each node-creating MPE
            for k = 1:length(obj.node.order)
                mpe = obj.mpe_by_name(obj.node.order(k).name);
                tt{k} = mpe.power_flow_node_types(obj, dm, obj.state.order(k).idx);
            end

            %% concatenate into a single node type vector
            ntv = vertcat(tt{:});

%             %% create node type struct
%             if nargout > 1
%                 %% define constants
%                 [PQ, PV, REF, NONE] = idx_bus;
% 
%                 ref = find(ntv == REF);     %% reference node indices
%                 pv  = find(ntv == PV );     %% PV node indices
%                 pq  = find(ntv == PQ );     %% PQ node indices
%                 nts = struct('ref', ref, 'pv', pv, 'pq', pq, ...
%                     'nref', length(ref), 'npv', length(pv), 'npq', length(pq));
%             end
        end

        function add_pf_constraints(obj, nm, om, dm, mpopt)
            %% system constraints
            obj.add_pf_system_constraints(om, dm, mpopt);
            
%             %% each element adds its PF constraints
%             for mpe = obj.mpe_list
%                 mpe{1}.add_pf_constraints(nm, om, dm, mpopt);
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
            for mpe = obj.mpe_list
                mpe{1}.add_opf_vars(nm, om, dm, mpopt);
            end
            
            %% legacy user-defined variables
            obj.add_opf_legacy_user_vars(om, dm, mpopt);
        end

        function add_opf_constraints(obj, nm, om, dm, mpopt)
            %% system constraints
            obj.add_opf_system_constraints(om, dm, mpopt);
            
            %% each element adds its OPF constraints
            for mpe = obj.mpe_list
                mpe{1}.add_opf_constraints(nm, om, dm, mpopt);
            end
        end

        function add_opf_costs(obj, nm, om, dm, mpopt)
            %% system costs
            obj.add_opf_system_costs(om, dm, mpopt);
            
            %% each element adds its OPF costs
            for mpe = obj.mpe_list
                mpe{1}.add_opf_costs(nm, om, dm, mpopt);
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
            %% create (read-only) copies of individual fields for convenience
            mpc = dm.mpc;
            [baseMVA, bus, gen, branch, gencost, Au, lbu, ubu, mpopt, ...
                N, fparm, H, Cw, z0, zl, zu, userfcn] = opf_args(mpc, mpopt);

            %% get some more problem dimensions
            if isfield(mpc, 'A')
                nlin = size(mpc.A, 1);  %% number of linear user constraints
            else
                nlin = 0;
            end
            if isfield(mpc, 'N')
                nw = size(mpc.N, 1);    %% number of general cost vars, w
            else
                nw = 0;
            end

            %% get number of user vars, check consistency
            nx = sum(cellfun(@(x)om.getN('var', x), obj.opf_legacy_user_var_names()));
            if nlin
                nz = size(mpc.A, 2) - nx; %% number of user z variables
                if nz < 0
                    error('mpe_network/add_opf_legacy_user_vars: user supplied A matrix must have at least %d columns.', nx);
                end
            else
                nz = 0;               %% number of user z variables
                if nw                 %% still need to check number of columns of N
                    if size(mpc.N, 2) ~= nx;
                        error('mpe_network/add_opf_legacy_user_vars: user supplied N matrix must have %d columns.', nx);
                    end
                end
            end

            %% save data
            om.userdata.user_vars = obj.opf_legacy_user_var_names();
            om.userdata.nlin = nlin;
            om.userdata.nw = nw;
            om.userdata.A = Au;
            om.userdata.l = lbu;
            om.userdata.u = ubu;
            om.userdata.N = N;
            om.userdata.fparm = fparm;
            om.userdata.H = H;
            om.userdata.Cw = Cw;

            %% add any user-defined vars
            if nz > 0
                om.add_var('z', nz, z0, zl, zu);
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
