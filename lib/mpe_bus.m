classdef mpe_bus < mp_element

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'bus';
%     end

    methods
        %% constructor
        function obj = mpe_bus()
            obj@mp_element();
            obj.name = 'bus';
            obj.np = 0;             %% this is a 0 port element
        end

        function obj = add_nodes(obj, nm, dm)
            dme = obj.data_model_element(dm);
            nm.add_node(obj.name, obj.nk);
        end

        %%-----  PF methods  -----
        function ntv = power_flow_node_types(obj, nm, dm, idx)
            %% define constants
            [PQ, PV, REF, NONE] = idx_bus;

            nb = obj.nk;            %% number of buses
            [ref, pv, pq] = bustypes(dm.mpc.bus, dm.mpc.gen);
            ntv = zeros(nb, 1);
            ntv(ref) = REF;
            ntv(pv)  = PV;
            ntv(pq)  = PQ;
        end
    end     %% methods
end         %% classdef
