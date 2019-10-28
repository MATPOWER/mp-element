classdef mp_aggregate < mp_element & mp_modeler
%MP_AGGREGATE

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        element_classes = {};   %% classes for individual element types
                                %% filled in by subclass init()
        mpe_list = {};          %% cell array of mp_element objects
        mpe_idx  = struct();    %% key = element name, val = index into mpe_list
        nv = 0;                 %% total number of v variables
        node = [];
        zvar = [];
        var = [];
    end
    
    methods
        %% constructor
        function obj = mp_aggregate(varargin)
%             fprintf('--> mp_aggregate (%d args)\n', nargin);
            obj@mp_element(varargin{:});
            if isempty(obj.node)    %% skip if constructed from existing object
                obj.init_set_types();   %% should be called in mp_modeler
                                        %% constructor, if not for:
                                        %% https://savannah.gnu.org/bugs/?52614
            end
            obj.name = 'aggregate';
            obj.mpc_field = '';
            obj.np = 0;     %% unknown number of ports at this point, init to 0
            obj.nk = 1;
            obj.nz = 0;     %% unknown number of z vars at this point, init to 0
%             fprintf('<-- mp_aggregate (%d args)\n', nargin);
        end

        function obj = create_model(obj, mpc)
            %% create element objects for each class with data
            for c = obj.element_classes
                mpe = c{1}();       %% element constructor
                if mpe.count(mpc)
                    obj.mpe_list{end+1} = mpe;
                    obj.mpe_idx.(mpe.name) = length(obj.mpe_list);
                    obj.np = obj.np + mpe.np * mpe.nk;  %% number of ports
                    obj.nz = obj.nz + mpe.nz * mpe.nk;  %% number of z vars
                end
            end
            
            %% create nodes
            obj.add_nodes(obj, mpc);
            
            %% add variables for each element type
            for mpe = obj.mpe_list
                mpe{1}.add_vars(obj, mpc);
            end

            %% build params
            for mpe = obj.mpe_list
                mpe{1}.build_params(obj, mpc);
            end
        end

        function mpe = mpe_by_name(obj, name)
            mpe = obj.mpe_list{obj.mpe_idx.(name)};
        end

        function obj = add_nodes(obj, asm, mpc)
            %% give each element opportunity to add nodes
            for mpe = obj.mpe_list
                mpe{1}.add_nodes(obj, mpc);
            end
            
            %% add voltage variables for each node
            obj.add_vvars(mpc);
        end

        function obj = add_vvars(obj, mpc)
            for k = 1:length(obj.node.order)
                mpe = obj.mpe_by_name(obj.node.order(k).name);
                mpe.add_vvars(obj, mpc);
            end
            obj.nv = obj.getN('var');
        end

        function obj = def_set_types(obj)
            obj.set_types = struct(...
                    'node', 'node', ...
                    'zvar', 'non-voltage state variable', ...
                    'var', 'variable' ...
                );
        end

        function obj = init_set_types(obj)
            %% call parent to create base data structures for each type
            init_set_types@mp_modeler(obj);

            %% finish initializing data structures for each type
            obj.node.data = struct( ...
                'idx2ID', struct(), ...
                'ID2idx', struct() );

            obj.zvar.data = struct( ...
                'idx2ID', struct(), ...
                'ID2idx', struct() );

            obj.var.data = struct( ...
                'v0', struct(), ...
                'vl', struct(), ...
                'vu', struct(), ...
                'vt', struct() );
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
                obj.node.data.v0 = subsasgn(obj.node.data.idx2ID, sc, idx2ID);
                obj.node.data.v0 = subsasgn(obj.node.data.ID2idx, sc, ID2idx);
            end
        end

        function add_zvar(obj, name, idx, varargin)
            %   obj.add_zvar(name, N, IDs)
            %   obj.add_zvar(name, N)
            %   obj.add_zvar(name, idx_list, N, IDs)
            %   obj.add_zvar(name, idx_list, N)
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
                    error('mp_aggregate/add_zvar: length of IDs vector (%d) must equal N (%d) ', length(idx2ID), N);
                end
            end
            if isempty(idx2ID)
                idx2ID = [1:N]';
            end

            %% create reverse mapping
            ID2idx = sparse(idx2ID, ones(N, 1), 1:N, max(idx2ID), 1);

            %% add the named zvar set
            obj.add_named_set('zvar', name, idx, N);
            
            %% add type-specific data for zvars (idx2ID, ID2idx)
            if isempty(idx)
                obj.zvar.data.idx2ID.(name) = idx2ID;
                obj.zvar.data.ID2idx.(name) = ID2idx;
            else
                %% calls to substruct() are relatively expensive, so we
                %% pre-build the struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                obj.zvar.data.v0 = subsasgn(obj.zvar.data.idx2ID, sc, idx2ID);
                obj.zvar.data.v0 = subsasgn(obj.zvar.data.ID2idx, sc, ID2idx);
            end
        end

        function add_var(obj, name, idx, varargin)
            %   obj.add_var(name, N, v0, vl, vu, vt)
            %   obj.add_var(name, N, v0, vl, vu)
            %   obj.add_var(name, N, v0, vl)
            %   obj.add_var(name, N, v0)
            %   obj.add_var(name, N)
            %   obj.add_var(name, idx_list, N, v0, vl, vu, vt)
            %   obj.add_var(name, idx_list, N, v0, vl, vu)
            %   obj.add_var(name, idx_list, N, v0, vl)
            %   obj.add_var(name, idx_list, N, v0)
            %   obj.add_var(name, idx_list, N)
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
            obj.add_named_set('var', name, idx, N);
            
            %% add type-specific data for var (v0, vl, vu, vt)
            if isempty(idx)
                obj.var.data.v0.(name) = v0;    %% initial value
                obj.var.data.vl.(name) = vl;    %% lower bound
                obj.var.data.vu.(name) = vu;    %% upper bound
                obj.var.data.vt.(name) = vt;    %% variable type
            else
                %% calls to substruct() are relatively expensive, so we
                %% pre-build the struct for addressing cell array fields
                %% sc = substruct('.', name, '{}', idx);
                sc = struct('type', {'.', '{}'}, 'subs', {name, idx});  %% cell array field
                obj.var.data.v0 = subsasgn(obj.var.data.v0, sc, v0);    %% initial value
                obj.var.data.vl = subsasgn(obj.var.data.vl, sc, vl);    %% lower bound
                obj.var.data.vu = subsasgn(obj.var.data.vu, sc, vu);    %% upper bound
                obj.var.data.vt = subsasgn(obj.var.data.vt, sc, vt);    %% variable type
            end
        end
    end     %% methods
end         %% classdef
