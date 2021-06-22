classdef dme_load_mpc2 < dme_load & dm_format_mpc2
%DME_LOAD_MPC2  MATPOWER data model load table for MATPOWER case format v2

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
        %% cannot be removed until we have a way to export dme_load to mpc
        function dm = parameterized(obj, dm, dmb, dmt, lam)
            dm = parameterized@dme_load(obj, dm, dmb, dmt, lam);    %% call parent

            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            cols = [PD QD];
            dm.mpc.bus(:, cols) = dmb.mpc.bus(:, cols) + ...
                lam * (dmt.mpc.bus(:, cols) - dmb.mpc.bus(:, cols));
        end
    end     %% methods
end         %% classdef
