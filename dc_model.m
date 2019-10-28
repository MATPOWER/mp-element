classdef dc_model < mp_model
%DC_MODEL  MATPOWER Model base class for DC models.
%   Each concrete MATPOWER Element class must inherit, at least indirectly,
%   from both MP_ELEMENT and MP_MODEL.
%
%   Subclass of MP_MODEL.
%   MP_MODEL provides propoerties and methods related to the specific
%   model and formulation (e.g. DC version, AC polar power version, etc.)
%
%   DC_MODEL defines:
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
            vtypes = {'v'};
        end
        function vtypes = model_zvars(obj)
            vtypes = {'z'};
        end

        function [B, K, p] = get_params(obj)
            np = obj.nk * obj.np;
            nz = obj.nk * obj.nz;
            if isempty(obj.B)
                B = sparse(np, np);
            else
                B = obj.B;
            end
            if isempty(obj.K)
                K = sparse(np, nz);
            else
                K = obj.K;
            end
            if isempty(obj.p)
                p = zeros(np, 1);
            else
                p = obj.p;
            end
        end

        function P = port_inj_power(obj, x, sysx, idx)
            % sys x : 1 = system x, 0 = class aggregate x
            
            [B, K, p] = obj.get_params();
            [v, z] = obj.x2vz(x, sysx);
            
            if sysx
                Ct = obj.C.';
                Dt = obj.D.';
                P = B*Ct*v + K*Dt*z + p;
            else
                P = B*v + K*z + p;
            end
        end
    end     %% methods
end         %% classdef
