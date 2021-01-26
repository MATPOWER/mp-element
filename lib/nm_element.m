classdef nm_element < handle
%NM_ELEMENT  Abstract base class for MATPOWER network model elements
%   NME = NM_ELEMENT()
%
%   Each concrete sub-class must also inherit from a sub-class of MP_FORM.
%
%   Properties
%       name : name of element type (constant across forumations)
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
%       count() - returns the number of elements of this type in dm, sets nme.nk
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
%       opf_add_vars()
%       opf_add_constraints()
%       opf_add_costs()

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        name = 'nm_element';
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
        function dme = data_model_element(obj, dm, name)
            if nargin < 3
                name = obj.name;
            end
            dme = dm.elm_by_name(name);
        end

        function nk = count(obj, dm)
            nk = dm.online(obj.name);
            obj.nk = nk;    %% update the count stored internally
        end

        function obj = add_nodes(obj, nm, dm)
        end

        function obj = add_states(obj, nm, dm)
        end

        function obj = add_vvars(obj, nm, dm, idx)
        end

        function obj = add_zvars(obj, nm, dm, idx)
        end

        function obj = build_params(obj, nm, dm)
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
            if isa(obj, 'mp_form')
                fprintf('FORMULATION NAME       : %s\n', obj.form_name());
                fprintf('FORMULATION TAG        : %s\n', obj.form_tag());
                fprintf('FORMULATION CLASS      : %s\n', obj.find_form_class());
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
%         function obj = pf_add_vars(obj, mm, nm, dm, mpopt)
%         end
% 
%         function obj = pf_add_constraints(obj, mm, nm, dm, mpopt)
%         end

        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
        end


        %%-----  CPF methods  -----
%         function cpf_add_vars(obj, mm, nm, dm, mpopt)
%         end
% 
%         function cpf_add_constraints(obj, mm, nm, dm, mpopt)
%         end


        %%-----  OPF methods  -----
        function obj = opf_add_vars(obj, mm, nm, dm, mpopt)
        end

        function obj = opf_add_constraints(obj, mm, nm, dm, mpopt)
        end

        function obj = opf_add_costs(obj, mm, nm, dm, mpopt)
        end

        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
        end
    end     %% methods
end         %% classdef
