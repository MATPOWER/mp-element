classdef dme_bus_nld_mpc2_node_test < dme_bus_mpc2
%DME_BUS_NLD_MPC2_NODE_TEST  MATPOWER data model bus table for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        bidx    %% indices into bus matrix (all rows) for "non-load buses"
    end     %% properties

    methods
        %% constructor
        function obj = dme_bus_nld_mpc2_node_test()
            obj@dme_bus_mpc2();     %% call parent constructor
            obj.name = 'bus_nld';
        end

        function nr = count(obj, dm)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS] = idx_bus;

            if isfield(dm.mpc, obj.table)
                tab = dm.mpc.(obj.table);
                obj.bidx = find(~tab(:, PD) & ~tab(:, QD) & ~tab(:, GS) & ~tab(:, BS));
                nr = length(obj.bidx);
            else
                nr = 0;
            end
            obj.nr = nr;
        end

        function tab = get_table(obj, dm)
            tab = dm.mpc.(obj.table)(obj.bidx, :);
        end

        function [gbus, ig] = gbus_vector(obj, gen_dme)
            %% buses of online gens
            ig = find(gen_dme.bus_type == 1);
            gbus = obj.i2on(gen_dme.bus(gen_dme.on(gen_dme.bus_type == 1)));
        end

        function midx = dme_idx2mpc_idx(obj, didx)
            if nargin > 1
                midx = obj.bidx(obj.on(didx));
            else
                midx = obj.bidx(obj.on);
            end
        end
    end     %% methods
end         %% classdef
