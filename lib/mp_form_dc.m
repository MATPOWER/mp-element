classdef mp_form_dc < mp_form
%MP_FORM_DC  MATPOWER Model base class for DC models.
%   Each concrete Network Model Element class must inherit, at least
%   indirectly, from both NM_ELEMENT and MP_FORM.
%
%   Subclass of MP_FORM.
%   MP_FORM provides properties and methods related to the specific
%   model and formulation (e.g. DC version, AC polar power version, etc.)
%
%   MP_FORM_DC defines:
%       linear active power injection = B theta + K z + p
%
%   Properties
%       (model parameters)
%       params - cell array of model parameter field names
%       B - np*nk x nn matrix
%       K - np*nk x nz matrix
%       p - np*nk x 1 matrix
%
%   Methods
%       model_name() - returns string w/name of model/formulation ('DC model')
%       model_tag() - returns string w/short label for model/formulation ('dc')
%       model_params() - cell array of names of model parameters
%                        {'B', 'K', 'p'}

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        %% model parameters
        B = [];
        K = [];
        p = [];
        param_ncols = struct('B', 2, 'K', 3, 'p', 1);
            %% num of columns for each parameter, where 1 = 1, 2 = np, 3 = nz
    end

    methods
        function name = model_name(obj)
            name = 'DC model';
        end
        function tag = model_tag(obj)
            tag = 'dc';
        end
        function params = model_params(obj)
           params = {'B', 'K', 'p'};
        end
        function vtypes = model_vvars(obj)
            vtypes = {'va'};
        end
        function vtypes = model_zvars(obj)
            vtypes = {'z'};
        end

        function P = port_inj_power(obj, x, sysx, idx)
            % sys x : 1 = system x, 0 = class aggregate x
            
            if nargin < 4
                idx = [];
            end
            
            [B, K, p] = obj.get_params(idx);
            [v, z] = obj.x2vz(x, sysx);
            
            if isempty(z)
                P = B*v + p;
            else
                P = B*v + K*z + p;
            end
%             if sysx
%                 Ct = obj.C';
%                 Dt = obj.D';
%                 if nargin < 4       %% all ports
%                     P = B*Ct*v + K*Dt*z + p;
%                 else                %% selected ports
%                     P = B(idx, :)*Ct*v + K(idx, :)*Dt*z + p(idx);
%                 end
%             else
%                 if nargin < 4       %% all ports
%                     P = B*v + K*z + p;
%                 else                %% selected ports
%                     P = B(idx, :)*v + K(idx, :)*z + p(idx);
%                 end
%             end
        end
    end     %% methods
end         %% classdef
