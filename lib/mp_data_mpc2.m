classdef mp_data_mpc2 < mp_data
%MP_DATA_MPC2  Implementation of MATPOWER data model for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpc         %% MATPOWER case struct
        user_mods   %% struct for storing input data for user variables,
                    %% constraints and costs
    end     %% properties

    methods
        function obj = build(obj, mpc, dmc)
            if nargin < 3
                dmc = [];
            end
            obj.mpc = loadcase(mpc);
            obj.baseMVA = obj.mpc.baseMVA;
            build@mp_data(obj, mpc, dmc);
        end

        function ntv = node_type_vector(obj, isref, ispv, ispq)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE] = idx_bus;

            %% package up bus type vector
            ntv = isref * REF + ispv * PV + ispq * PQ;
        end

        function ref = node_type_ref(obj, ntv)
            %% define constants
            [PQ, PV, REF, NONE] = idx_bus;

            ref = find(ntv == REF);
        end

        function pv = node_type_pv(obj, ntv)
            %% define constants
            [PQ, PV, REF, NONE] = idx_bus;

            pv = find(ntv == PV);
        end

        function pq = node_type_pq(obj, ntv)
            %% define constants
            [PQ, PV, REF, NONE] = idx_bus;

            pq = find(ntv == PQ);
        end

%         function display(obj)
%             display@mp_data(obj);
%             mpc = obj.mpc
%         end

        function print_soln(obj, fname)
        end

        function save_soln(obj, fname)
        end
    end     %% methods
end         %% classdef
