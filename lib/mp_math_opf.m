classdef mp_math_opf < mp_math
%MP_MATH_OPF  MATPOWER mathematical model for optimal power flow (OPF) problem.
%   ?
%
%   MP_MATH_OPF ... optimal power flow ...
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
%     end     %% properties

    methods
        function obj = build(obj, nm, dm, mpopt)
            nm.opf_add_vars(obj, nm, dm, mpopt);
            nm.opf_add_constraints(obj, nm, dm, mpopt);
            nm.opf_add_costs(obj, nm, dm, mpopt);
        end

        function dm = data_model_update(obj, nm, dm, mpopt)
            nm.opf_data_model_update(obj, nm, dm, mpopt);
        end

        function nm = network_model_x_soln(obj, nm)
            %% convert solved state from math model to network model soln
            [nm.soln.v, nm.soln.z, nm.soln.x] = ...
                nm.opf_convert_x(obj.soln.x, obj.get_userdata('aux_data'));
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            opt = nm.opf_solve_opts(obj, dm, mpopt);
        end
    end     %% methods
end         %% classdef
