classdef nme_shunt_ac < nme_shunt% & mp_form_ac

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
        function obj = build_params(obj, nm, dm)
            build_params@nme_shunt(obj, nm, dm);    %% call parent

            dme = obj.data_model_element(dm);
            Ysh = dme.Gs + 1j * dme.Bs;             %% shunt admittances
            nsh = obj.nk;
            obj.Y = sparse(1:nsh, 1:nsh, Ysh, nsh, nsh);
        end

        %%-----  CPF methods  -----
        function obj = cpf_data_model_update(obj, mm, nm, dm, mpopt)
            ad = mm.get_userdata('aux_data');
            dme = obj.data_model_element(dm);
            dm = dme.parameterized(dm, ad.dmb, ad.dmt, mm.soln.x(end));
        end
    end     %% methods
end         %% classdef
