classdef mp_aggregate < mp_element & mp_idx_manager% & mp_model
%MP_AGGREGATE Abstract class, explicitly a subclass of mp_element and
%             mp_idx_manager and implicitly assumed to be subclasses of
%             mp_model as well

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
        function obj = mp_aggregate()
            obj@mp_element();
            obj.name = 'aggregate';
            obj.mpc_field = '';
            obj.np = 0;     %% unknown number of ports at this point, init to 0
            obj.nk = 1;
            obj.nz = 0;     %% unknown number of z_ vars at this point, init to 0
        end

        function obj = create_model(obj, mpc, mpopt)
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
                mpc.sys_wide_zip_loads = mpopt.exp.sys_wide_zip_loads;
            end
            %%-----  end of HACK  -----

            %% create element objects for each class with data
            i = 0;
            for c = obj.element_classes
                mpe = c{1}();       %% element constructor
                if mpe.count(mpc)
                    i = i + 1;
                    obj.mpe_list{i} = mpe;
                    obj.mpe_map.(mpe.name) = i;
                    obj.np = obj.np + mpe.np * mpe.nk;  %% number of ports
                    obj.nz = obj.nz + mpe.nz * mpe.nk;  %% number of z_ vars
                end
            end
            
            if obj.np ~= 0      %% skip for empty model
                %% create nodes and node voltage state variables
                obj.add_nodes(obj, mpc);
            
                %% create non-voltage states and corresponding state variables
                obj.add_states(obj, mpc);
            
                %% build params
                obj.build_params(obj, mpc);
            end
        end

        function mpe = mpe_by_name(obj, name)
            mpe = obj.mpe_list{obj.mpe_map.(name)};
        end

        function obj = add_nodes(obj, asm, mpc)
            %% each element adds its nodes
            for mpe = obj.mpe_list
                mpe{1}.add_nodes(obj, mpc);
            end
            
            %% add voltage variables for each node
            obj.add_vvars(obj, mpc);
        end

        function obj = add_states(obj, asm, mpc)
            %% each element adds its states
            for mpe = obj.mpe_list
                mpe{1}.add_states(obj, mpc);
            end
            
            %% add state variables for each node
            obj.add_zvars(obj, mpc);
        end

        function obj = build_params(obj, asm, mpc)
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
                mpe.build_params(obj, mpc);
                C = horzcat(C, mpe.C);
                D = horzcat(D, mpe.D);
                pk = pk + mpe.np * mpe.nk;          %% increment port counter
                zk = zk + mpe.nz * mpe.nk;          %% increment z-var counter
                obj.mpe_port_map(k, 2) = pk - 1;    %% ending port index
                obj.mpe_z_map(k, 2)    = zk - 1;    %% ending z-var index
            end
            obj.C = { horzcat(C{:}) };
            obj.D = { horzcat(D{:}) };
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

        function obj = add_vvars(obj, asm, mpc, idx)
            for k = 1:length(obj.node.order)
                mpe = obj.mpe_by_name(obj.node.order(k).name);
                mpe.add_vvars(obj, mpc, obj.state.order(k).idx);
            end
            for vtype = obj.model_vvars
                obj.nv = obj.nv + obj.getN(vtype{1});
            end
        end

        function obj = add_zvars(obj, asm, mpc, idx)
            for k = 1:length(obj.state.order)
                mpe = obj.mpe_by_name(obj.state.order(k).name);
                mpe.add_zvars(obj, mpc, obj.state.order(k).idx);
            end
        end

        %%-----  mp_idx_manager methods  -----
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
                    error('mp_aggregate/add_node: length of IDs vector (%d) must equal N (%d) ', length(idx2ID), N);
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
                    error('mp_aggregate/add_state: length of IDs vector (%d) must equal N (%d) ', length(idx2ID), N);
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
                if nargs >= 2
                    vl = args{2};
                    if nargs >= 3
                        vu = args{3};
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

        %%-----  OPF methods  -----
        function add_opf_vars(obj, asm, om, mpc, mpopt)
            vars = horzcat(obj.model_vvars(), obj.model_zvars());
            for vtype = vars
                st = obj.(vtype{1});
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
                mpe{1}.add_opf_vars(asm, om, mpc, mpopt);
            end
        end

        function add_opf_constraints(obj, asm, om, mpc, mpopt)
            %% system constraints
            obj.add_opf_system_constraints(om, mpc, mpopt);
            
            %% each element adds its OPF constraints
            for mpe = obj.mpe_list
                mpe{1}.add_opf_constraints(asm, om, mpc, mpopt);
            end
        end

        function add_opf_costs(obj, asm, om, mpc, mpopt)
            %% system costs
            obj.add_opf_system_costs(om, mpc, mpopt);
            
            %% each element adds its OPF costs
            for mpe = obj.mpe_list
                mpe{1}.add_opf_costs(asm, om, mpc, mpopt);
            end
        end

        function add_opf_system_constraints(obj, om, mpc, mpopt)
            %% can be overridden to add additional system constraints
            %% node balance constraints
            obj.add_opf_node_balance_constraints(om);
        end

        function add_opf_system_costs(obj, om, mpc, mpopt)
            %% can be overridden to add additional system costs
        end
    end     %% methods
end         %% classdef
