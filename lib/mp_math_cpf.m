classdef mp_math_cpf < mp_math_pf
%MP_MATH_CPF  MATPOWER math model for continuation power flow (CPF) problem.
%   ?
%
%   MP_MATH_CPF ... continuation power flow ...
%
%   Properties
%       ? - ?
%
%   Methods
%       ?

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
        function obj = add_aux_data(obj, nm, dm, mpopt)
            %% create aux_data struct
            obj.aux_data = obj.cpf_aux_data(nm, dm, mpopt);
        end

        function ad = cpf_aux_data(obj, nm, dm, mpopt)
            dmt = dm.userdata.target;
            nmt = nm.userdata.target;

            ad  = obj.pf_aux_data(nm, dm, mpopt);
            adt = obj.pf_aux_data(nmt, dmt, mpopt);

            ad = nm.cpf_check_xfer(nmt, ad, adt);
            ad.nmt = nmt;
            ad.dmb = dm.copy(); %% save copy of unmodified base case data model
            ad.dmt = dmt;
            ad.adt = adt;
        end

        function add_vars(obj, nm, dm, mpopt)
            add_vars@mp_math_pf(obj, nm, dm, mpopt);
            %% required for using nmt to compute nm state from mm state
            obj.aux_data.adt.var_map = obj.aux_data.var_map;
            obj.add_var('lambda', 1, 0);
        end

        function dm = data_model_update(obj, nm, dm, mpopt)
            nm.cpf_data_model_update(obj, nm, dm, mpopt);
        end

        function nm = network_model_x_soln(obj, nm)
            %% convert solved state from math model to network model soln
            [nm.soln.v, nm.soln.z, nm.soln.x] = ...
                nm.cpf_convert_x(obj.soln.x, obj.aux_data);
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            opt = nm.cpf_solve_opts(obj, dm, mpopt);
        end
    end     %% methods
end         %% classdef
