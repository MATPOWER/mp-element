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
        bus     %% bus index vector (all gens)
        Pg0     %% initial active power (p.u.) for gens that are on
        Qg0     %% initial reactive power (p.u.) for gens that are on
        Vg      %% generator voltage setpoint
        Pmin    %% active power lower bound (p.u.) for gens that are on
        Pmax    %% active power upper bound (p.u.) for gens that are on
        Qmin    %% reactive power lower bound (p.u.) for gens that are on
        Qmax    %% reactive power upper bound (p.u.) for gens that are on
        pcost   %% active power cost parameters for gens that are on
        qcost   %% reactive power cost parameters for gens that are on
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
