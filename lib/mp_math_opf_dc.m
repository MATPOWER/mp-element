classdef mp_math_opf_dc < mp_math_opf
%MP_MATH_OPF_DC  MATPOWER mathematical model for DC optimal power flow (OPF) problem.
%   ?
%
%   MP_MATH_OPF_DC ... power flow ...
%
%   Properties
%       ? - ?
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        %% constructor
        function obj = mp_math_opf_dc()
            obj@mp_math_opf();
            obj.element_classes = { @mme_bus_opf_dc, @mme_gen_opf_dc, @mme_branch_opf_dc };
        end

        function tag = form_tag(obj)
            tag = 'dc';
        end

        function name = form_name(obj)
            name = 'DC';
        end

        function [vx, z, x] = opf_convert_x(obj, mmx, nm)
            nm_vars = obj.update_nm_vars(mmx, nm);

            %% convert (real) math model x to network model x
            vx = nm_vars.va;
            z  = nm_vars.z;
            if nargout < 2
                vx = [vx; z];
            elseif nargout > 2
                x = [vx; z];
            end
        end

        function names = legacy_user_var_names(obj)
            names = {'Va', 'Pg'};
        end

        function obj = add_node_balance_constraints(obj, nm, dm, mpopt)
            [B, K, p] = nm.get_params();

            %% power balance constraints
            C = nm.C;
            Amis = C * [B*C' K*nm.D'];
            bmis = -C * p;
            obj.add_lin_constraint('Pmis', Amis, bmis, bmis, ...
                                [nm.va.order; nm.z.order]);
        end

        function add_system_costs(obj, nm, dm, mpopt)
            %% legacy user-defined costs
            if isfield(dm.userdata, 'legacy_opf_user_mods')
                obj.add_legacy_user_costs(nm, dm, 1);
            end
        end

        function opt = solve_opts(obj, nm, dm, mpopt)
            opt = mpopt2qpopt(mpopt, obj.problem_type());

            switch opt.alg
                case {'MIPS', 'IPOPT'}
                    if mpopt.opf.start < 2      %% initialize interior point
                        opt.x0 = obj.interior_x0(obj, nm, dm);
                    end
                case 'OSQP'
                    opt.x0 = [];        %% disable provided starting point
            end
        end
    end     %% methods
end         %% classdef
