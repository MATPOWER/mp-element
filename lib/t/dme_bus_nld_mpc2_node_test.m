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
        bus_class   %% 1 = bus_nld, 2 = bus_ld
        bidx    %% indices into bus matrix (all rows) for this type of bus
    end     %% properties

    methods
        %% constructor
        function obj = dme_bus_nld_mpc2_node_test()
            obj@dme_bus_mpc2();     %% call parent constructor
            obj.bus_class = 1;
            obj.name = 'bus_nld';
        end

        function bidx = mpc_idx(obj, tab)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS] = idx_bus;

            bidx = find(~tab(:, PD) & ~tab(:, QD) & ~tab(:, GS) & ~tab(:, BS));
        end

        function nr = count(obj, dm)
            if isfield(dm.mpc, obj.table)
                tab = dm.mpc.(obj.table);
                obj.bidx = obj.mpc_idx(tab);
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
            ig = find(gen_dme.bus_type == obj.bus_class);
            gbus = obj.i2on(gen_dme.bus(gen_dme.on(gen_dme.bus_type == obj.bus_class)));
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
