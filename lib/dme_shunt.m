classdef dme_shunt < dm_element
%DME_SHUNT  Abstract base class for MATPOWER data model shunt table

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bus     %% bus index vector (all shunts)
        Gs      %% shunt conductance (p.u. active power demanded at
                %% V = 1.0 p.u.) for shunts that are on
        Bs      %% shunt susceptance (p.u. reactive power injected at
                %% V = 1.0 p.u.) for shunts that are on
    end     %% properties

    methods
        %% constructor
        function obj = dme_shunt()
            obj@dm_element();   %% call parent constructor
            obj.name = 'shunt';
            obj.table = 'bus';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus', 'gs', 'bs', 'p', 'q'});
        end
    end     %% methods
end         %% classdef
