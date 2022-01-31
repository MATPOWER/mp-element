classdef mp_math_opf_ac < mp_math_opf
%MP_MATH_OPF_AC  MATPOWER mathematical model for AC optimal power flow (OPF) problem.
%   ?
%
%   MP_MATH_OPF_AC ... power flow ...
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
        function add_system_costs(obj, nm, dm, mpopt)
            %% legacy user-defined costs
            if isfield(dm.userdata, 'legacy_opf_user_mods')
                obj.add_legacy_user_costs(nm, dm, 0);
            end
        end

        function add_legacy_user_constraints(obj, nm, dm, mpopt)
            %% call parent
            add_legacy_user_constraints@mp_math_opf(obj, nm, dm, mpopt);

            if ~isempty(dm.userdata.legacy_opf_user_mods)
                uc = dm.userdata.legacy_opf_user_mods.nlc;
                for k = 1:length(uc)
                    obj.add_nln_constraint(uc{k}{:});
                end
            end
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            opt = mpopt2nlpopt(mpopt, obj.problem_type());

            if mpopt.opf.start < 2
                %% initialize interior point
                opt.x0 = nm.opf_interior_x0(obj, nm, dm);
            end
        end
    end     %% methods
end         %% classdef
