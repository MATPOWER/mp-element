classdef dme_shunt < dm_element
%DME_SHUNT  Abstract base class for MATPOWER data model shunt table

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        %% constructor
        function obj = dme_shunt()
            obj@dm_element();   %% call parent constructor
            obj.name = 'shunt';
            obj.table = 'bus';
        end
    end     %% methods
end         %% classdef