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
            b2i = dm.elements.bus.ID2i;     %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            tab = obj.get_table(dm);
            obj.fbus = b2i(tab(:, F_BUS));
            obj.tbus = b2i(tab(:, T_BUS));
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elements.bus.status;    %% bus status

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

            branch = obj.get_table(dm);

            obj.R  = branch(obj.on, BR_R);
            obj.X  = branch(obj.on, BR_X);
            obj.B  = branch(obj.on, BR_B);
            obj.tap    = branch(obj.on, TAP);
            obj.shift  = branch(obj.on, SHIFT) * pi/180;
            obj.rate_a = branch(obj.on, RATE_A) / dm.baseMVA;
        end

        function obj = update(obj, dm, varargin)
            %% obj.update(dm, name1, val1, name2, val2, ...)
            %% obj.update(dm, idx, name1, val1, name2, val2, ...)

            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
            baseMVA = dm.baseMVA;

            n = length(varargin);
            if rem(n, 2)    %% odd
                idx = obj.on(varargin{1});
                s = 2;      %% starting arg index
            else            %% even
                idx = obj.on;
                s = 1;      %% starting arg index
            end
            for k = s:2:n-1
                val = varargin{k+1};
                switch varargin{k}
                    case 'Sf'
                        dm.mpc.branch(idx, PF) = real(val) * baseMVA;
                        dm.mpc.branch(idx, QF) = imag(val) * baseMVA;
                    case 'St'
                        dm.mpc.branch(idx, PT) = real(val) * baseMVA;
                        dm.mpc.branch(idx, QT) = imag(val) * baseMVA;
                    case 'Pf'
                        dm.mpc.branch(idx, PF) = val * baseMVA;
                    case 'Pt'
                        dm.mpc.branch(idx, PT) = val * baseMVA;
                    case 'Qf'
                        dm.mpc.branch(idx, QF) = val * baseMVA;
                    case 'Qt'
                        dm.mpc.branch(idx, QT) = val * baseMVA;
                    case {'muSf', 'muPf'}
                        dm.mpc.branch(idx, MU_SF) = val / baseMVA;
                    case {'muSt', 'muPt'}
                        dm.mpc.branch(idx, MU_ST) = val / baseMVA;
                    case 'muAngmin'
                        dm.mpc.branch(idx, MU_ANGMIN) = val * pi/180;
                    case 'muAngmax'
                        dm.mpc.branch(idx, MU_ANGMAX) = val * pi/180;
                end
            end
        end
    end     %% methods
end         %% classdef
