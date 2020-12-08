classdef dme_bus_mpc2 < dme_bus & dm_format_mpc2
%DME_BUS  MATPOWER data model bus table for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties
    methods
        %% constructor
        function obj = dme_bus_mpc2()
            obj@dme_bus();      %% call parent constructor

            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;
            obj.id_col = BUS_I;
            obj.st_col = BUS_TYPE;
        end

        function status = get_status(obj, dm)
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;

            %% check that all buses have a valid BUS_TYPE
            tab = obj.get_table(dm);
            bt = tab(:, obj.st_col);
            err = find(~(bt == PQ | bt == PV | bt == REF | bt == NONE));
            if ~isempty(err)
                error('dme_bus_mpc2/get_status: bus %d has an invalid BUS_TYPE', err);
            end
            status = (bt ~= NONE);       %% bus status
            obj.status = status;
        end
    end     %% methods
end         %% classdef
