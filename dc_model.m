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
    end     %% methods
end         %% classdef
