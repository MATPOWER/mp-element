classdef mpe_branch < mp_element

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'branch';
%     end
    
    methods
        %% constructor
        function obj = mpe_branch()
            obj@mp_element();
            obj.name = 'branch';
            obj.np = 2;             %% this is a 2 port element
        end

        function obj = build_params(obj, nm, dm)
            dme = obj.data_model_element(dm);

            %% incidence matrices
            nidx = nm.get_node_idx('bus');  %% node indices for 'bus'
            fidx = nidx(dme.fbus);          %% branch "from" port node indices
            tidx = nidx(dme.tbus);          %% branch "to" port node indices
            obj.C = obj.incidence_matrix(nm.getN('node'), fidx, tidx);
            obj.D = obj.incidence_matrix(nm.getN('state'));
        end
    end     %% methods
end         %% classdef
