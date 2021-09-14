classdef mp_math < opt_model
%MP_MATH  MATPOWER mathematical model abstract base class.
%   ?
%
%   MP_MATH provides properties and methods related to the specific
%   problem specification being solved (e.g. power flow, continuation
%   power flow, optimal power flow, etc.) ...
%
%   Properties
%       ? - ?
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        aux_data    %% struct of auxiliary data relevant to the model,
                    %% e.g. can be passed to model constraint functions
    end

    methods
        function display(obj)
            fprintf('MATH MODEL CLASS : %s\n', class(obj));
            display@opt_model(obj)
        end
    end     %% methods
end         %% classdef
