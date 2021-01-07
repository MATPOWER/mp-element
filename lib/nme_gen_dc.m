classdef nme_gen_dc < nme_gen & mp_form_dc

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'gen';
%     end
    
    methods
        function obj = add_zvars(obj, nm, dm, idx)
            ng = obj.nk;
            dme = obj.data_model_element(dm);
            nm.add_var('z', 'Pg', ng, dme.Pg0, dme.Pmin, dme.Pmax);
        end

        function obj = build_params(obj, nm, dm)
            build_params@nme_gen(obj, nm, dm);      %% call parent
            ng = obj.nk;
            obj.K = -speye(ng);
        end

        %%-----  OPF methods  -----
        function opf_build_gen_cost_params(obj, dm)
            dme = obj.data_model_element(dm);
            obj.cost = dme.opf_build_gen_cost_params(dm, 1);
        end

        function obj = opf_add_constraints(obj, mm, nm, dm, mpopt)
            %% piecewise linear costs
            if obj.cost.pwl.n
                mm.add_lin_constraint('ycon', obj.cost.pwl.A, [], obj.cost.pwl.b, {'Pg', 'y'});
            end

            %% call parent
            opf_add_constraints@nme_gen(obj, mm, nm, dm, mpopt);
        end

        function obj = opf_data_model_update(obj, mm, nm, dm, mpopt)
            %% generator active power
            ss = nm.get_idx('state');
            Pg = nm.soln.z(ss.i1.gen:ss.iN.gen);

            %% shadow prices on generator limits
            vv = mm.get_idx();
            lambda = mm.soln.lambda;
            muPmax = lambda.upper(vv.i1.Pg:vv.iN.Pg);
            muPmin = lambda.lower(vv.i1.Pg:vv.iN.Pg);

            %% update in the data model
            dme = obj.data_model_element(dm);
            dme.update(dm, Pg, muPmin, muPmax);
        end
    end     %% methods
end         %% classdef
