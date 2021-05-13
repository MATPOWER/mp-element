classdef nme_bus_ld_acp_node_test < nme_bus_nld_acp_node_test

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
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
        function obj = nme_bus_ld_acp_node_test()
            obj@nme_bus_nld_acp_node_test();
            obj.name = 'bus_ld';
            obj.np = 1;             %% this is a 1 port element
        end

        function obj = build_params(obj, nm, dm)
            %% incidence matrices
            nidx = nm.get_node_idx('bus_ld');   %% node indices for 'bus_ld'
            idx = nidx([1:obj.nk]');            %% node indices for loads
            obj.C = obj.incidence_matrix(nm.getN('node'), idx);
            obj.D = obj.incidence_matrix(nm.getN('state'));

            dme = obj.data_model_element(dm);

            %% bus shunts
            Ysh = dme.Gs + 1j * dme.Bs;             %% shunt admittances
            obj.Y = sparse(1:obj.nk, 1:obj.nk, Ysh, obj.nk, obj.nk);

            %% constant power loads
            obj.s = dme.Pd + 1j * dme.Qd;           %% complex power demand
        end
    end     %% methods
end         %% classdef
