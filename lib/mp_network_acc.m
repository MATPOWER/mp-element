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


        %%-----  OPF methods  -----
        function [vx_, z_, x_] = opf_convert_x(obj, mmx, ad)
            %% convert (real) math model x to (complex) network model x_
            nv_ = obj.nv / 2;       %% number of voltage vars (sysx=1)
            nz_ = obj.nz;           %% number of state vars
            vr = mmx(1:nv_, :);     b = nv_;
            vi = mmx(b+1:b+nv_, :); b = b + nv_;
            zr = mmx(b+1:b+nz_, :); b = b + nz_;
            zi = mmx(b+1:b+nz_, :);
            vx_ = vr+1j*vi;
            z_  = zr+1j*zi;
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
