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
        function obj = mp_math_opf_dc()
            obj@mp_math_opf();
            obj.element_classes = { @mme_opf_gen_dc, @mme_opf_branch_dc };
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
    end     %% methods
end         %% classdef
