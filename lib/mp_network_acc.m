classdef mp_network_acc < mp_network_ac% & mp_form_acc

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        vr = [];
        vi = [];
    end
    
    methods
        function obj = mp_network_acc()
            obj@mp_network_ac();
            obj.element_classes = ...
                { @nme_bus_acc, @nme_gen_acc, @nme_load_acc, ...
                    @nme_branch_acc, @nme_shunt_acc };

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  INIT_SET_TYPES() is called explicitly as needed
            %%              (if obj.node is empty) in BUILD() and DISPLAY(),
            %%              after object construction, but before object use.
        end

        function obj = def_set_types(obj)
            def_set_types@mp_network_ac(obj);   %% call parent first
            obj.set_types.vr = 'REAL VOLTAGE VARS (vr)';
            obj.set_types.vi = 'IMAG VOLTAGE VARS (vi)';
        end

        function va = initial_voltage_angle(obj, idx)
            vr = obj.params_var('vr');  %% inital value
            vi = obj.params_var('vi');  %% inital value
            if nargin < 2 || isempty(idx)
                va = angle(vr + 1j * vi);
            else
                va = angle(vr(idx) + 1j * vi(idx));
            end
        end


        %%-----  PF methods  -----
        function [vx_, z_, x_] = pf_convert_x(obj, mmx, ad, only_v)
            %% x = obj.pf_convert(mmx, ad)
            %% [v, z] = obj.pf_convert(mmx, ad)
            %% [v, z, x] = obj.pf_convert(mmx, ad)
            %% ... = obj.pf_convert(mmx, ad, only_v)

            %% update v_, z_ from mmx
            nm_vars = obj.update_vars(mmx, ad);
            vx_ = nm_vars.vr + 1j * nm_vars.vi;
            z_ = nm_vars.zr + 1j * nm_vars.zi;

            %% update z, if requested
            if nargin < 4 || ~only_v
                z_ = obj.pf_update_z(vx_, z_, ad);
            end

            if nargout < 2
                vx_ = [vx_; z_];
            elseif nargout > 2
                x_ = [vx_; z_];
            end
        end


        %%-----  OPF methods  -----
        function [vx_, z_, x_] = opf_convert_x(obj, mmx, ad)
            nm_vars = obj.update_vars(mmx, ad);

            %% convert (real) math model x to (complex) network model x_
            vx_ = nm_vars.vr + 1j * nm_vars.vi;
            z_  = nm_vars.zr + 1j * nm_vars.zi;
            if nargout < 2
                vx_ = [vx_; z_];
            elseif nargout > 2
                x_ = [vx_; z_];
            end
        end

        function names = opf_legacy_user_var_names(obj)
            names = {'Vr', 'Vi', 'Pg', 'Qg'};
        end
    end     %% methods
end         %% classdef
