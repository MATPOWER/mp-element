classdef mp_network < nm_element & mpe_container & mp_idx_manager% & mp_form
%MP_NETWORK  Abstract base class for MATPOWER network model
%   Explicitly a subclass of NM_ELEMENT, MP_IDX_MANAGER and MPE_CONTAINER,
%   and implicitly assumed to be a subclass of MP_FORM as well.

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        nv = 0;                 %% total number of (real) v variables
        node = [];
        port = [];
        state = [];
    end
    
    methods
        %% constructor
        function obj = mp_network()
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
            %%              (if obj.node is empty) in BUILD() and DISPLAY(),
            %%              after object construction, but before object use.
        end

        function obj = build(obj, dm)
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
            for k = 1:length(obj.elm_list)
                nme = obj.elm_list{k};
                nme.build_params(obj, dm);
                C = horzcat(C, {nme.C});
                D = horzcat(D, {nme.D});

                %% add ports (to track indexing of all network ports)
                if nme.np > 1
                    obj.init_indexed_name('port', nme.name, {nme.np});
                    for p = 1:nme.np
                        obj.add_port(nme.name, {p}, nme.nk);
                    end
                elseif nme.np
                    obj.add_port(nme.name, nme.nk);
                end
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
                    'state', 'STATES', ...
                    'port', 'PORTS' ...
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

        function add_port(obj, name, idx, N)
            %   obj.add_port(name, N)
            %   obj.add_port(name, idx_list, N)
            if ~iscell(idx)
                N = idx;
                idx = {};
            end

            %% add the named state set
            obj.add_named_set('port', name, idx, N);
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

        function s = set_type_idx_map(obj, set_type, idxs, dm, group_by_name)
            %%  s = obj.set_type_idx_map(set_type, idxs)
            %%  s = obj.set_type_idx_map(set_type, idxs, dm)
            %%  s = obj.set_type_idx_map(set_type, idxs, dm, group_by_name)
            %%  output struct s has fields:
            %%      name - name of named set
            %%      idx - optional cell array of indices of named set
            %%      i - index of element in named set
            %%      e - external index (i.e. corresp. row in data model)
            %%      ID - external ID (i.e. corresp. element ID in data model)
            %%      j - (only for group_by_name == 1), corresponding index of
            %%          set type (i.e. node, port or state)
            s = set_type_idx_map@mp_idx_manager(obj, set_type, idxs);
            n = length(idxs(:));
            if nargin > 3
                for k = 1:n
                    dme = dm.elm_by_name(s(k).name);
                    s(k).e = dme.on(s(k).i);        %% external index
                    if isempty(dme.ID)
                        s(k).ID = s(k).e;
                    else
                        s(k).ID = dme.ID(s(k).e);   %% external ID
                    end
                end

                %% consolidate indices into vectors for each unique
                %% name/idx pair, if requested
                if nargin > 4 && group_by_name
                    %% extract fields
                    [name, idx, i, e, ID] = deal(cell(size(idxs)));
                    [name{:}] = deal(s.name);
                    [idx{:}]  = deal(s.idx);
                    [i{:}]  = deal(s.i);    i = cell2mat(i);
                    [e{:}]  = deal(s.e);    e = cell2mat(e);
                    [ID{:}] = deal(s.ID);   ID = cell2mat(ID);

                    %% find unique name/idx
                    name_idx = cellfun(@join_name_idx, name, idx, ...
                        'UniformOutput', 0);
                    [c, ia, ic] = unique(name_idx);

                    %% recreate struct, grouped by name/idx
                    c0 = cell(size(c));
                    s = struct('name', name(ia), 'idx', idx(ia), ...
                        'i', c0, 'e', c0, 'ID', c0, 'j', c0);
                    for k = 1:length(ia)
                        s(k).i  = i(ic == k);
                        s(k).e  = e(ic == k);
                        s(k).ID = ID(ic == k);
                        s(k).j = idxs(ic == k);
                    end
                end
            end
        end

        function label = set_type_label(obj, set_type, idxs, dm)
            %%  label = obj.set_type_label(set_type, idxs)
            %%  label = obj.set_type_label(set_type, idxs, dm)
            label = cell(size(idxs));
            if nargin > 3
                s = obj.set_type_idx_map(set_type, idxs, dm);
                [ID{1:length(idxs(:))}] = deal(s.ID);
            else
                s = obj.set_type_idx_map(set_type, idxs);
                [ID{1:length(idxs(:))}] = deal(s.i);
            end
            for k = 1:length(idxs(:))
                if isempty(s(k).idx)
                    label{k} = sprintf('%s %d', s(k).name, ID{k});
                else
                    if length(s(k).idx) <= 1
                        idxstr = sprintf('%d', s(k).idx{1});
                    else
                        idxstr = [sprintf('%d', s(k).idx{1}) sprintf(',%d', s(k).idx{2:end})];
                    end
                    label{k} = sprintf('%s(%s) %d', s(k).name, idxstr, ID{k});
                end
            end
            if isscalar(idxs)               %% return scalar
                label = label{1};
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

        function varargout = get_node_idx(obj, name)
            %% [i1 iN] = obj.get_node_idx(name)
            %% nidx = obj.get_node_idx(name), where
            %%      nidx = [i1:iN]' or
            %%      nidx = {[i1(1):iN(1)]', ..., [i1(n):iN(n)]'}
            [varargout{1:nargout}] = obj.get_set_type_idx('node', name);
        end

        function varargout = get_port_idx(obj, name)
            %% [i1 iN] = obj.get_port_idx(name)
            %% pidx = obj.get_port_idx(name), where
            %%      pidx = [i1:iN]' or
            %%      pidx = {[i1(1):iN(1)]', ..., [i1(n):iN(n)]'}
            [varargout{1:nargout}] = obj.get_set_type_idx('port', name);
        end

        function varargout = get_state_idx(obj, name)
            %% [i1 iN] = obj.get_state_idx(name)
            %% sidx = obj.get_state_idx(name), where
            %%      sidx = [i1:iN]' or
            %%      sidx = {[i1(1):iN(1)]', ..., [i1(n):iN(n)]'}
            [varargout{1:nargout}] = obj.get_set_type_idx('state', name);
        end

        function [i1, iN] = get_set_type_idx(obj, node_state, name)
            %% private method
            %% [i1 iN] = obj.get_set_type_idx(node_state, name)
            %% nidx = obj.get_set_type_idx(node_state, name), where
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
            %%           ntv = obj.node_types(nm, dm)
            %% [ref, pv, pq] = obj.node_types(nm, dm)
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
        function obj = pf_add_constraints(obj, mm, nm, dm, mpopt)
            %% system constraints
            obj.pf_add_system_constraints(mm, dm, mpopt);
            
%             %% each element adds its PF constraints
%             for nme = obj.elm_list
%                 nme{1}.pf_add_constraints(mm, nm, dm, mpopt);
%             end
        end

        function obj = pf_add_system_constraints(obj, mm, dm, mpopt)
            %% can be overridden to add additional system constraints

            %% node balance constraints
            obj.pf_add_node_balance_constraints(mm, dm, mpopt);
        end

        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            %% each element updates its data model
            for nme = obj.elm_list
                nme{1}.pf_data_model_update(mm, nm, dm, mpopt);
            end
        end

        opt = pf_solve_opts(obj, mm, dm, mpopt)


        %%-----  OPF methods  -----
        function obj = opf_soln(obj, mm)
        end
        
        function obj = opf_add_vars(obj, mm, nm, dm, mpopt)
            %% add network voltage and non-voltage state variables
            vars = horzcat(obj.model_vvars(), obj.model_zvars());
            for vtype = vars
                st = obj.(vtype{1});    %% set type
                for k = 1:st.NS
                    name = st.order(k).name;
                    if isempty(st.order(k).idx)
                        d = st.data;
                        mm.add_var(name, st.idx.N.(name), d.v0.(name), d.vl.(name), d.vu.(name), d.vt.(name));
                    else
                        error('mp_network/opf_add_vars: handling of indexed sets not implmented here (yet)');
                    end
                end
            end
            
            %% each element adds its OPF variables
            for nme = obj.elm_list
                nme{1}.opf_add_vars(mm, nm, dm, mpopt);
            end
            
            %% legacy user-defined variables
            obj.opf_add_legacy_user_vars(mm, dm, mpopt);
        end

        function obj = opf_add_constraints(obj, mm, nm, dm, mpopt)
            %% system constraints
            obj.opf_add_system_constraints(mm, dm, mpopt);
            
            %% each element adds its OPF constraints
            for nme = obj.elm_list
                nme{1}.opf_add_constraints(mm, nm, dm, mpopt);
            end
        end

        function obj = opf_add_costs(obj, mm, nm, dm, mpopt)
            %% system costs
            obj.opf_add_system_costs(mm, dm, mpopt);
            
            %% each element adds its OPF costs
            for nme = obj.elm_list
                nme{1}.opf_add_costs(mm, nm, dm, mpopt);
            end
        end

        function opf_add_system_constraints(obj, mm, dm, mpopt)
            %% can be overridden to add additional system constraints

            %% node balance constraints
            obj.opf_add_node_balance_constraints(mm);

            %% legacy user-defined constraints
            obj.opf_add_legacy_user_constraints(mm, dm, mpopt);
        end

        function opf_add_system_costs(obj, mm, dm, mpopt)
        end

        function opf_add_legacy_user_vars(obj, mm, dm, mpopt)
            z = dm.user_mods.z;

            %% save data
            mm.userdata.user_vars = obj.opf_legacy_user_var_names();

            %% add any user-defined vars
            if z.nz > 0
                mm.add_var('z', z.nz, z.z0, z.zl, z.zu);
                mm.userdata.user_vars{end+1} = 'z';
            end
        end

        function opf_add_legacy_user_constraints(obj, mm, dm, mpopt)
            lin = dm.user_mods.lin;

            %% user-defined linear constraints
            if lin.nlin
                uv = mm.get_userdata('user_vars');
                mm.add_lin_constraint('usr', lin.A, lin.l, lin.u, uv);
            end
        end

        function opf_add_legacy_user_costs(obj, mm, dm, dc)
            user_cost = dm.user_mods.cost;
            if user_cost.nw
                uv = mm.get_userdata('user_vars');
                mm.add_legacy_cost('usr', user_cost, uv);
            end

            %% implement legacy user costs using quadratic or general non-linear costs
            cp = mm.params_legacy_cost();   %% construct/fetch the parameters
            [N, H, Cw, rh, m] = deal(cp.N, cp.H, cp.Cw, cp.rh, cp.mm);
            [nw, nx] = size(N);
            if nw
                if any(cp.dd ~= 1) || any(cp.kk)    %% not simple quadratic form
                    if dc                           %% (includes "dead zone" or
                        if any(cp.dd ~= 1)          %%  quadratic "penalty")
                            error('mp_network/opf_add_legacy_user_costs: DC OPF can only handle legacy user-defined costs with d = 1');
                        end
                        if any(cp.kk)
                            error('mp_network/opf_add_legacy_user_costs: DC OPF can only handle legacy user-defined costs with no "dead zone", i.e. k = 0');
                        end
                    else
                        %% use general nonlinear cost to implement legacy user cost
                        legacy_cost_fcn = @(x)opf_legacy_user_cost_fcn(x, cp);
                        mm.add_nln_cost('usr', 1, legacy_cost_fcn);
                    end
                else                                %% simple quadratic form
                    %% use a quadratic cost to implement legacy user cost
                    %% f = 1/2 * w'*H*w + Cw'*w, where w = diag(m)*(N*x - rh)
                    %% Let: MN = diag(m)*N
                    %%      MR = M * rh
                    %%      HMR  = H  * MR;
                    %%      HtMR = H' * MR;
                    %%  =>   w = MN*x - MR
                    %% f = 1/2 * (MN*x - MR)'*H*(MN*x - MR) + Cw'*(MN*x - MR)
                    %%   = 1/2 * x'*MN'*H*MN*x +
                    %%          (Cw'*MN - 1/2 * MR'*(H+H')*MN)*x +
                    %%          1/2 * MR'*H*MR - Cw'*MR
                    %%   = 1/2 * x'*Q*w + c'*x + k

                    M    = sparse(1:nw, 1:nw, m, nw, nw);
                    MN   = M * N;
                    MR   = M * rh;
                    HMR  = H  * MR;
                    HtMR = H' * MR;
                    Q = MN' * H * MN;
                    c = full(MN' * (Cw - 1/2*(HMR+HtMR)));
                    k = (1/2 * HtMR - Cw)' * MR;
                    mm.add_quad_cost('usr', Q, c, k);
                end
            end
        end

        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
            %% each element updates its data model
            for nme = obj.elm_list
                nme{1}.opf_data_model_update(mm, nm, dm, mpopt);
            end
        end
    end     %% methods
end         %% classdef

function name_idx = join_name_idx(name, idx)
    if isempty(idx)
        name_idx = name;
    else
        name_idx = [name sprintf('_%d', idx{:})];
    end
end
