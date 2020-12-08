classdef dme_branch < dm_element
%DME_BRANCH  Abstract base class for MATPOWER data model branch table

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        fbus        %% bus index vector for "from" port (port 1)
        tbus        %% bus index vector for "to" port (port 2)
    end     %% properties

    methods
        %% constructor
        function obj = dme_branch()
            obj@dm_element();   %% call parent constructor
            obj.name = 'branch';
            obj.table = 'branch';
        end
    end     %% methods
end         %% classdef
