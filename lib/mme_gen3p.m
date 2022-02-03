classdef mme_gen3p < mm_element

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'gen3p';
%     end
    
    methods
        %% constructor
        function obj = mme_gen3p()
            obj@mm_element();
            obj.name = 'gen3p';
        end

        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            dme = obj.data_model_element(dm);
            nme = obj.network_model_element(nm);

            ss = nm.get_idx('state');

            for p = 1:nme.nz
                %% generator active power
                sg = nm.soln.z(ss.i1.gen3p(p):ss.iN.gen3p(p)) * dm.base_kva;
                pg = real(sg);

                %% update in the data model
                dme.tab.(sprintf('pg%d', p))(dme.on) = real(sg);
                dme.tab.(sprintf('pf%d', p))(dme.on) = cos(angle(sg));
            end
        end
    end     %% methods
end         %% classdef
