classdef mp_math_pf < mp_math
%MP_MATH_PF  MATPOWER mathematical model for power flow (PF) problem.
%   ?
%
%   MP_MATH_PF ... power flow ...
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
        function obj = build(obj, nm, dm, mpopt)
            build@mp_math(obj, nm, dm, mpopt);  %% call parent
            obj.add_aux_data(nm, dm, mpopt);
            obj.add_vars(nm, dm, mpopt);
            obj.add_constraints(nm, dm, mpopt);
        end

        function obj = add_aux_data(obj, nm, dm, mpopt)
            %% create aux_data struct
            obj.aux_data = obj.pf_aux_data(nm, dm, mpopt);
        end

        function obj = add_system_vars(obj, nm, dm, mpopt)
            obj.add_pf_system_vars(nm, dm, mpopt);
        end

        function dm = data_model_update(obj, nm, dm, mpopt)
            nm.pf_data_model_update(obj, nm, dm, mpopt);
        end

        function nm = network_model_x_soln(obj, nm)
            [nm.soln.v, nm.soln.z, nm.soln.x] = ...
                obj.pf_convert_x(obj.soln.x, nm, obj.aux_data);
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            %% for AC power flow only, must override for mp_math_pf_dc
            switch mpopt.pf.alg
                case 'DEFAULT'
                    opt = mpopt2nleqopt(mpopt, obj.problem_type(), 'DEFAULT');
                case {'NR', 'NR-SP', 'NR-SC', 'NR-SH', 'NR-IP', 'NR-IC', 'NR-IH'}
                    opt = mpopt2nleqopt(mpopt, obj.problem_type(), 'NEWTON');
                case {'FDXB', 'FDBX'}
                    opt = mpopt2nleqopt(mpopt, obj.problem_type(), 'FD');
                    opt.fd_opt.jac_approx_fcn = @()nm.pf_fd_jac_approx(obj, dm, mpopt);
                    opt.fd_opt.labels = {'P', 'Q'};
                case 'FSOLVE'
                    opt = mpopt2nleqopt(mpopt, obj.problem_type(), 'FSOLVE');
                case 'GS'
                    opt = mpopt2nleqopt(mpopt, obj.problem_type(), 'GS');
                    opt.gs_opt.x_update_fcn = ...
                        @(x, f)obj.gs_x_update(x, f, nm, dm, mpopt);
                case 'ZG'
                    opt = mpopt2nleqopt(mpopt, obj.problem_type(), 'ZG');
                    zg_x_update = @(x, f)obj.zg_x_update(x, f, nm, dm, mpopt);
                    opt.core_sp = struct(...
                        'alg',              'ZG', ...
                        'name',             'Implicit Z-bus Gauss', ...
                        'default_max_it',   1000, ...
                        'need_jac',         0, ...
                        'update_fcn',       zg_x_update  );
                otherwise
                    error('mp_math_pf/solve_opts: invalid value for MPOPT.PF.ALG (%s)', mpopt.pf.alg);
            end
            opt.verbose = mpopt.verbose;
        end
    end     %% methods
end         %% classdef
