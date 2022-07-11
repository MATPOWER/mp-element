classdef (Abstract) mp_math_opf_acc < mp_math_opf_ac
%MP_MATH_OPF_ACC  MATPOWER mathematical model for AC optimal power flow (OPF) problem.
%   ?
%
%   MP_MATH_OPF_ACC ... power flow ...
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
        function obj = mp_math_opf_acc()
            obj@mp_math_opf_ac();
            obj.element_classes = { @mme_bus_opf_acc, @mme_gen_opf_ac, ...
                @mme_load_pf_ac, @mme_branch_opf_acc, @mme_shunt_pf_ac };
        end

        function [vx_, z_, x_] = convert_x_m2n(obj, mmx, nm)
            nm_vars = obj.update_nm_vars(mmx, nm);

            %% convert (real) math model x to (complex) network model x_
            vx_ = nm_vars.vr + 1j * nm_vars.vi;
            z_  = nm_vars.zr + 1j * nm_vars.zi;
            if nargout < 2
                vx_ = [vx_; z_];
            elseif nargout > 2
                x_ = [vx_; z_];
            end
        end

        function varef1 = interior_va(obj, nm, dm)
            %% return scalar va equal to angle of first reference node
            ad = obj.aux_data;
            ref1 = ad.ref(1);
            varef1 = angle(ad.vr(ref1) + 1j * ad.vi(ref1));
        end
    end     %% methods
end         %% classdef
