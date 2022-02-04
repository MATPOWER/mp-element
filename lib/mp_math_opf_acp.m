classdef mp_math_opf_acp < mp_math_opf_ac
%MP_MATH_OPF_ACP  MATPOWER mathematical model for AC optimal power flow (OPF) problem.
%   ?
%
%   MP_MATH_OPF_ACP ... power flow ...
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
        function obj = mp_math_opf_acp()
            obj@mp_math_opf_ac();
            obj.element_classes = { @mme_bus_opf_acp, @mme_gen_opf_ac, ...
                @mme_branch_opf_acp, @mme_bus3p_opf_acp, @mme_gen3p_opf, ...
                @mme_line3p_opf, @mme_xfmr3p_opf, @mme_buslink_opf_acp };
        end

        function [vx_, z_, x_] = opf_convert_x(obj, mmx, nm)
            nm_vars = obj.update_nm_vars(mmx, nm);

            %% convert (real) math model x to (complex) network model x_
            vx_ = nm_vars.vm .* exp(1j * nm_vars.va);
            z_  = nm_vars.zr + 1j * nm_vars.zi;
            if nargout < 2
                vx_ = [vx_; z_];
            elseif nargout > 2
                x_ = [vx_; z_];
            end
        end

        function names = legacy_user_var_names(obj)
            names = {'Va', 'Vm', 'Pg', 'Qg'};
        end
    end     %% methods
end         %% classdef
