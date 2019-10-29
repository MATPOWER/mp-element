classdef mp_element < handle
%MP_ELEMENT  MATPOWER Element abstract base class.
%   MPE = MP_ELEMENT(FORMULATION, MPC)
%   MPE = MP_ELEMENT(MPE0)
%   MPE = MP_ELEMENT(S)
%
%   Properties
%       np : number of ports per element
%       nz : number of non-voltage state variables per element
%       nk : number of elements
%       C : cell array of sparse element-node incidence matrices,
%           where C{j}(i,k) is 1 if port j of element k is connected to node i
%       D : cell array of sparse incidence matrices for Z variables,
%           where D{j} is for Z variable j
%       
%
%   Methods
%       constructor - should be called explicitly without args at top of
%                     sub-class constructors
%           mpe = mp_element(mpc)
%           mpe = mp_element(other object)
%       init() - called automatically by constructor with same args, should
%           be called explicitly by sub-class constructor if called with args
%           mpe.init(...)
%       add_nodes() - add node, voltage variables, and node to element ID mappings
%       add_states() - add state, state variables
%           
%       build_params() - build model parameters from data model
%           params = mpe.build_params(asm, mpc)
%       count() - returns the number of elements of this type in mpc, sets mpe.nk
%           TorF = mpe.count(mpc)
%
%       S = mpe.port_inj_power(x, idx)
%       [S, dS] = mpe.port_inj_power(x, idx)
%       I = mpe.port_inj_current(x, idx)
%       [I, dI] = mpe.port_inj_current(x, idx)

%       [nn, n2b, b2n] = mpe.nodes()
%           nn  : number of nodes
%           n2b : (internal) node number to (external) bus number mapping
%           b2n : (external) bus number to (internal) node number mapping
%
%       [params, varsets] = mpe.power_balance(form)
%           form : polar vs. rectangular
%           params : struct of parameters needed to compute power balance
%               depending on form
%               e.g. polar form g_s(x) = S_bus(V) + C*x + d + g_u(x)
%                    where S_bus(V) = diag(V) * conj(sum_over_k(Ybus_k)*V)
%               Ybus_k : Ybus term from this mpe
%               C_k : C term from this mpe
%               d_k : d term from this mpe
%               g_k : g_u term from this mpe
%               dg_k : Jacobian of g_u term from this mpe
%               d2G_k : Hessian of g_u term from this mpe
%           varsets : variable sets corresponding to C_k and g_k, etc.
%
%       newvars = mpe.variables()
%           newvars
%               name
%               N
%               type : state vs. control?

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        name = 'mp_element';
        mpc_field = '';     %% field checked for presence of this element type
        np = 0;             %% number of ports per element
        nk = 0;             %% number of elements of this type loaded
        nz = 0;             %% number of non-voltage states per element
        C = {};             %% cell array of sparse element-node incidence
                            %% matrices, where C{j}(i,k) is 1 if port j of
                            %% element k is connected to node i
        D = {};             %% cell array of sparse incidence matrices for
                            %% Z variables, where D{j} is for Z variable j
    end
    
    methods
        %% constructor
        %% subclass constructor should be identical (update parent class name)
        function obj = mp_element(varargin)
%             fprintf('--> mp_element (%d args)\n', nargin);
            if nargin > 0
                s = varargin{1};
                if isa(s, 'mp_element')         %% copy another object
                    props = fieldnames(s);
                    for k = 1:length(props)
                        obj.(props{k}) = s.(props{k});
                    end
                elseif ~isstruct(s) || ~isfield(s, 'bus')   %% check MPC
                    error('@mp_element/mp_element: input must be an ''mp_element'' object or MPC');
                end
            end
%             fprintf('<-- mp_element (%d args)\n', nargin);
        end

        function nk = count(obj, mpc)
            if isfield(mpc, obj.mpc_field) && ~isempty(mpc.(obj.mpc_field))
                nk = size(mpc.(obj.mpc_field), 1);
            end
            obj.nk = nk;        %% update the count stored internally
        end

        function obj = add_nodes(obj, asm, mpc)
        end

        function obj = add_states(obj, asm, mpc)
        end

        function obj = add_vvars(obj, asm, mpc)
        end

        function obj = add_zvars(obj, asm, mpc)
        end

        function obj = build_params(obj, asm, mpc)
        end

        function [v, z] = x2vz(obj, x, sysx);
            % sys x : 1 = system x, 0 = class aggregate x
            if sysx
                nv = size(obj.C{1}, 1);
%                 nz = size(obj.D{1}, 1);     % doesn't work when D = {}
            else
                nv = obj.nk * obj.np;
%                 nz = obj.nk * obj.nz;
            end
            v = x(1:nv);
            z = x(nv+1:end);
%             z = x(nv+1:nv+nz);
        end
    end     %% methods
end         %% classdef
