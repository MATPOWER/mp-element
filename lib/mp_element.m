classdef mp_element < handle
%MP_ELEMENT  MATPOWER Element abstract base class.
%   MPE = MP_ELEMENT()
%
%   Each concrete sub-class must also inherit from a sub-class of MP_MODEL.
%
%   Properties
%       name : name of element type (constant across forumations)
%       dm_table : name of table in data model to be checked for this
%           element type
%       np : number of ports per element
%       nz : number of non-voltage state variables per element
%       nk : number of elements
%       C : cell array of sparse element-node incidence matrices,
%           where C{j}(i,k) is 1 if port j of element k is connected to node i
%       D : cell array of sparse incidence matrices for Z variables,
%           where D{j}(i,k) is 1 if j-th Z variable for element k corresponds
%           to element i of system Z
%
%   Methods
%       count() - returns the number of elements of this type in mpc, sets mpe.nk
%       get_nv_()
%       x2vz()
%       incidence_matrix()
%       display()
%
%   Abstract Methods
%       add_nodes()
%       add_states()
%       add_vvars()
%       add_zvars()
%       build_params() - build model parameters from data model
%       add_opf_vars()
%       add_opf_constraints()
%       add_opf_costs()

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        name = 'mp_element';
        dm_table = '';      %% name of table in dm to be checked for presence
                            %% of this element type
        np = 0;             %% number of ports per element
        nz = 0;             %% number of non-voltage states per element (possibly complex)
        nk = 0;             %% number of elements of this type loaded
        C = [];             %% stacked element-node incidence matrices,
                            %% where C(i,kk) is 1 if port j of element k is
                            %% connected to node i, and kk = k + (j-1)*np
        D = [];             %% stacked sparse incidence matrices for
                            %% Z variables, where D(i,kk) is 1 if z-variable j
                            %% of element k is the i-th system z-variable
                            %% and kk = k + (j-1)*nz
        soln                %% struct for storing solved states, quantities
    end
    
    methods
        function nk = count(obj, mpc)
            if isfield(mpc, obj.dm_table) && ~isempty(mpc.(obj.dm_table))
                nk = size(mpc.(obj.dm_table), 1);
                obj.nk = nk;    %% update the count stored internally
            else
                nk = 0;
            end
        end

        function obj = add_nodes(obj, nm, mpc)
        end

        function obj = add_states(obj, nm, mpc)
        end

        function obj = add_vvars(obj, nm, mpc, idx)
        end

        function obj = add_zvars(obj, nm, mpc, idx)
        end

        function obj = build_params(obj, nm, mpc)
        end

        function nv_ = get_nv_(obj, sysx);
            % sysx : 1 = system x_, 0 = element class x_

            %% get sizes
            if sysx
                nv_ = size(obj.C, 1);
%                 nz_ = size(obj.D, 1);
            else
                nv_ = obj.nk * obj.np;
%                 nz_ = obj.nk * obj.nz;
            end
        end

        function [v_, z_, vi_] = x2vz(obj, x_, sysx, idx);
            % sysx : 1 = system x_, 0 = element class x_
            % if x_ is a matrix, each output will have the same number of
            % columns, each column considered a separate instance of the vectors

            %% split x_
            nv_ = obj.get_nv_(sysx);
            v_ = x_(1:nv_, :);
            z_ = x_(nv_+1:end, :);

            %% set full port voltages and states for element class
            if sysx         %% system x_ is provided, convert to x_ for ports
                v_ = obj.C' * v_;   %% full port voltages for element class
                z_ = obj.D' * z_;   %% full states for element class
            end

            %% port voltages for selected ports
            if nargout > 2
                if isempty(idx)
                    vi_ = v_;
                else
                    vi_ = v_(idx, :);
                end
            end
        end

        function CD = incidence_matrix(obj, m, varargin)
            %% obj.incidence_matrix(m, idx1, idx2, ...)
            n = length(varargin);   %% number of ports/z-vars
            if n == 1
                CD = sparse(varargin{1}, 1:obj.nk, 1, m, obj.nk);
            elseif n > 1
                blocks = cell(1, n);
                for i = 1:n
                    blocks{i} = sparse(varargin{i}, 1:obj.nk, 1, m, obj.nk);
                end
                CD = horzcat(blocks{:});
            else
                CD = sparse(m, 0);
            end
        end

%         function A = getA(obj, tr)
%             if nargin < 1
%                 C = obj.C';
%                 D = obj.D';
%             else
%                 C = obj.C;
%                 D = obj.D;
%             end
%             [mC, nC] = size(C);
%             [mD, nD] = size(D);
%             A = [ C sparse(mC, nD); sparse(mD, nC) D ];
%         end
%         
%         function Ap = getAprime(obj, tr)
%             if nargin < 1
%                 C = obj.C';
%                 D = obj.D';
%             else
%                 C = obj.C;
%                 D = obj.D;
%             end
%             [mC, nC] = size(C);
%             [mD, nD] = size(D);
%             Ap = [   C sparse(mC, nC+2*nD);
%                     sparse(mC, nC) C sparse(mC, 2*nD);
%                     sparse(mD, 2*nC) D sparse(mD, nD);
%                     sparse(mD, 2*nC+nD) D ];
%         end
        
        function display(obj)
%             if have_feature('octave')
%                 struct(obj)
%             else
%                 display@handle(obj)
%             end
            fprintf('ELEMENT NAME           : %s\n', obj.name);
            fprintf('ELEMENT CLASS          : %s\n', class(obj));
            fprintf('# OF ELEMENTS          : %d\n', obj.nk);
            fprintf('# OF PORTS/ELEM        : %d\n', obj.np);
            fprintf('# OF NON-V STATES/ELEM : %d\n', obj.nz);
            if isa(obj, 'mp_model')
                fprintf('MODEL NAME             : %s\n', obj.model_name());
                fprintf('MODEL TAG              : %s\n', obj.model_tag());
                fprintf('MODEL CLASS            : %s\n', obj.find_model_class());
                fprintf('MODEL PARAMETERS');
                model_params = obj.model_params();
                for j = 1:length(model_params)
                    pn = model_params{j};   %% parameter name
                    if j == 1
                        fmt = '%6s : ';
                    else
                        fmt = '%22s : ';
                    end
                    if isempty(obj.(pn))
                        fprintf([fmt '-\n'], pn);
                    else
                        [m, n] = size(obj.(pn));
                        if ~full(any(any(obj.(pn))))
                            s = '(all zeros)';
                        else
                            s = '';
                        end
                        fprintf([fmt '%d x %-7d%s\n'], pn, m, n, s);
                    end
                end
            end
        end

        %%-----  PF methods  -----
%         function add_pf_vars(obj, nm, om, mpc, mpopt)
%         end
% 
%         function add_pf_constraints(obj, nm, om, ad, mpc, mpopt)
%         end


        %%-----  OPF methods  -----
        function add_opf_vars(obj, nm, om, mpc, mpopt)
        end

        function add_opf_constraints(obj, nm, om, mpc, mpopt)
        end

        function add_opf_costs(obj, nm, om, mpc, mpopt)
        end
    end     %% methods
end         %% classdef
