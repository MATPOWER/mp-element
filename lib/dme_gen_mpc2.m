classdef dme_gen_mpc2 < dme_gen & dm_format_mpc2
%DME_GEN_MPC2  MATPOWER data model gen table for MATPOWER case format v2

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
        function obj = dme_gen_mpc2()
            obj@dme_gen();      %% call parent constructor

            %% define constants
            [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS] = idx_gen;
            obj.st_col = GEN_STATUS;
        end

        function obj = initialize(obj, dm)
            initialize@dme_gen(obj, dm);    %% call parent

            %% define constants
            [GEN_BUS] = idx_gen;

            %% get bus mapping info
            b2i = dm.elm_by_name('bus').ID2i;   %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            tab = obj.get_table(dm);
            obj.bus = b2i(tab(:, GEN_BUS));
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elm_by_name('bus').status;  %% bus status

            %% update status of gens at isolated/offline buses
            obj.status = obj.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@dme_gen(obj, dm);
        end
    end     %% methods
end         %% classdef
