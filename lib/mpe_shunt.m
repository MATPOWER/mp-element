classdef mpe_shunt < mp_element

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        %% constructor
        function obj = mpe_shunt()
            obj@mp_element();
            obj.name = 'shunt';
            obj.np = 1;             %% this is a 1 port element
        end

        function obj = build_params(obj, nm, dm)
            dme = obj.data_model_element(dm);

            %% incidence matrices
            IDs = dme.busID(dme.on);                %% bus IDs
            nidx = nm.node.data.ID2idx.bus(IDs);    %% node indexes
            obj.C = obj.incidence_matrix(nm.getN('node'), nidx);
            obj.D = obj.incidence_matrix(nm.getN('state'));
        end
    end     %% methods
end         %% classdef
