classdef nme_gen3p < nm_element % & mp_form_ac

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        %% constructor
        function obj = nme_gen3p()
            obj@nm_element();
            obj.name = 'gen3p';
            obj.np = 3;             %% this is a 3 port element
            obj.nz = 3;
        end

        function obj = add_zvars(obj, nm, dm, idx)
            p = idx{1};
            ng = obj.nk;
            dme = obj.data_model_element(dm);
            pg_start = dme.(sprintf('pg%d_start', p));
            qg_start = dme.(sprintf('qg%d_start', p));

            if p == 1
                nm.init_indexed_name('zr', 'Pg', {obj.nz});
                nm.init_indexed_name('zi', 'Qg', {obj.nz});
            end
            nm.add_var('zr', 'Pg', {p}, ng, pg_start, 0, Inf);
            nm.add_var('zi', 'Qg', {p}, ng, qg_start, -Inf, Inf);
        end

        function obj = build_params(obj, nm, dm)
            build_params@nm_element(obj, nm, dm);   %% call parent
            obj.N = -speye(obj.nk * obj.nz);
        end

        %%-----  PF methods  -----
        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            dme = obj.data_model_element(dm);
            ss = nm.get_idx('state');

            for p = 1:obj.nz
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
