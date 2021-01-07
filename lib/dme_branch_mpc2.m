classdef dme_branch_mpc2 < dme_branch & dm_format_mpc2
%DME_BRANCH_MPC2  MATPOWER data model branch table for MATPOWER case format v2

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
        function obj = dme_branch_mpc2()
            obj@dme_branch();   %% call parent constructor

            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                    TAP, SHIFT, BR_STATUS] = idx_brch;
            obj.st_col = BR_STATUS;
        end

        function obj = initialize(obj, dm)
            initialize@dme_branch(obj, dm);     %% call parent

            %% define named indices into data matrices
            [F_BUS, T_BUS] = idx_brch;

            %% get bus mapping info
            b2i = dm.elm_by_name('bus').ID2i;   %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            tab = obj.get_table(dm);
            obj.fbus = b2i(tab(:, F_BUS));
            obj.tbus = b2i(tab(:, T_BUS));
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elm_by_name('bus').status;  %% bus status

            %% update status of branches connected to isolated/offline buses
            obj.status = obj.status & bs(obj.fbus) & ...
                                      bs(obj.tbus);

            %% call parent to fill in on/off
            update_status@dme_branch(obj, dm);
        end

        function obj = build_params(obj, dm)
            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
            baseMVA = dm.mpc.baseMVA;

            branch = obj.get_table(dm);

            obj.R  = branch(obj.on, BR_R);
            obj.X  = branch(obj.on, BR_X);
            obj.B  = branch(obj.on, BR_B);
            obj.tap    = branch(obj.on, TAP);
            obj.shift  = branch(obj.on, SHIFT) * pi/180;
            obj.rate_a = branch(obj.on, RATE_A) / baseMVA;
        end

        function obj = update(obj, dm, Sf, St, muSf, muSt, muAngmin, muAngmax)
            %% obj.update(dm, Sf, St)
            %% obj.update(dm, Sf, St, muSf, muSt)
            %% obj.update(dm, Sf, St, muSf, muSt, muAngmin, muAngmax)

            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
            baseMVA = dm.mpc.baseMVA;
 
            dm.mpc.branch(obj.on, PF) = real(Sf) * baseMVA;
            dm.mpc.branch(obj.on, PT) = real(St) * baseMVA;
            dm.mpc.branch(obj.on, QF) = imag(Sf) * baseMVA;
            dm.mpc.branch(obj.on, QT) = imag(St) * baseMVA;

            if nargin > 4
                dm.mpc.branch(obj.on, MU_SF) = muSf / baseMVA;
                dm.mpc.branch(obj.on, MU_ST) = muSt / baseMVA;
                if nargin > 6
                    dm.mpc.branch(obj.on, MU_ANGMIN) = muAngmin * pi/180;
                    dm.mpc.branch(obj.on, MU_ANGMAX) = muAngmax * pi/180;
                end
            end
        end
    end     %% methods
end         %% classdef
