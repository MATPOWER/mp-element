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
        function [vx_, z_, x_] = opf_convert_x(obj, mmx, nm)
            nm_vars = nm.update_vars(mmx, obj.aux_data);

            %% convert (real) math model x to (complex) network model x_
            vx_ = nm_vars.vm .* exp(1j * nm_vars.va);
            z_  = nm_vars.zr + 1j * nm_vars.zi;
            if nargout < 2
                vx_ = [vx_; z_];
            elseif nargout > 2
                x_ = [vx_; z_];
            end
        end

        function names = opf_legacy_user_var_names(obj)
            names = {'Va', 'Vm', 'Pg', 'Qg'};
        end
    end     %% methods
end         %% classdef
