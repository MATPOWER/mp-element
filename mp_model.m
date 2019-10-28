classdef mp_model < handle
%MP_MODEL  MATPOWER Model abstract base class.
%   Each concrete MATPOWER Element class must inherit, at least indirectly,
%   from both MP_ELEMENT and MP_MODEL.
%
%   MP_MODEL provides propoerties and methods related to the specific
%   model and formulation (e.g. DC version, AC polar power version, etc.)
%
%   Properties
%       subclasses provide properties for model parameters
%
%   Methods
%       model_name() - returns string w/name of model/formulation
%       model_tag() - returns string w/short label for model/formulation
%       model_params() - cell array of names of model parameters

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         mp_model_field = '';
%     end
    
    methods
        function name = model_name(obj)
            error('model_name() method not implemented');
        end
        function tag = model_tag(obj)
            error('model_tag() method not implemented');
        end
        function params = model_params(obj)
            error('model_params() method not implemented');
        end
    end     %% methods
end         %% classdef
