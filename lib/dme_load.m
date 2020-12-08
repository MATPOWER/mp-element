classdef dme_load < dm_element
%DME_LOAD  Abstract base class for MATPOWER data model load table

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus         %% bus index vector
        Pd
        Qd
    end     %% properties

    methods
        %% constructor
        function obj = dme_load()
            obj@dm_element();   %% call parent constructor
            obj.name = 'load';
            obj.table = 'bus';
        end
    end     %% methods
end         %% classdef
