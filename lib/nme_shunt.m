classdef nme_shunt < nm_element

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
        function obj = nme_shunt()
            obj@nm_element();
            obj.name = 'shunt';
            obj.np = 1;             %% this is a 1 port element
        end

        function obj = build_params(obj, nm, dm)
            dme = obj.data_model_element(dm);

            %% incidence matrices
            nidxs = obj.node_indices(nm, dm, 'bus', 'bus');
            obj.C = obj.incidence_matrix(nm.getN('node'), nidxs{:});
            obj.D = obj.incidence_matrix(nm.getN('state'));
        end
    end     %% methods
end         %% classdef
