classdef dme_gen < dm_element
%DME_GEN  Abstract base class for MATPOWER data model gen table

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        busID
    end     %% properties

    methods
        %% constructor
        function obj = dme_gen()
            obj@dm_element();   %% call parent constructor
            obj.name = 'gen';
            obj.table = 'gen';
        end
    end     %% methods
end         %% classdef
