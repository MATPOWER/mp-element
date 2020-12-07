classdef dme_load_mpc2 < dme_load & dm_format_mpc2
%DME_LOAD_MPC2  MATPOWER data model load table for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        busidx
    end     %% properties

    methods
        function nr = count(obj, dm)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD] = idx_bus;

            tab = obj.get_table(dm);
            obj.busidx = find(tab(:, PD) | tab(:, QD));
            nr = length(obj.busidx);
            obj.nr = nr;
        end

%         function obj = update_status(obj, dm)
%             %% define constants
%             [GEN_BUS] = idx_gen;
%             
%             dm_bus = dm.elm_by_name('bus');
%             bs = dm_bus.status;     %% bus status
%             b2i = dm_bus.ID2i;      %% bus num to idx mapping
% 
%             %% update status of loads at isolated/offline buses
%             tab = obj.get_table(dm);
%             obj.status = obj.status & bs(b2i(tab(:, GEN_BUS)));
% 
%             %% call parent to fill in on/off
%             update_status@dme_load(obj);
%         end
    end     %% methods
end         %% classdef
