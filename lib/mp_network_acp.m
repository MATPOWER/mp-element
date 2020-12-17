classdef mp_network_acp < mp_network_ac & mp_form_acp

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        va = [];
        vm = [];
    end
    
    methods
        function obj = mp_network_acp()
            obj@mp_network_ac();
            obj.element_classes = ...
                { @nme_bus_acp, @nme_gen_acp, @nme_load_acp, ...
                    @nme_branch_acp, @nme_shunt_acp };

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
            obj.set_types.va = 'VOLTAGE ANG VARS (va)';
            obj.set_types.vm = 'VOLTAGE MAG VARS (vm)';
        end


        %%-----  OPF methods  -----
        function x_ = opf_convert_x(obj, mmx, ad)
            %% convert (real) math model x to (complex) network model x_
            nv_ = obj.nv / 2;       %% number of voltage vars (sysx=1)
            nz_ = obj.nz;           %% number of state vars
            va = mmx(1:nv_, :);     b = nv_;
            vm = mmx(b+1:b+nv_, :); b = b + nv_;
            zr = mmx(b+1:b+nz_, :); b = b + nz_;
            zi = mmx(b+1:b+nz_, :);
            x_ = [vm .* exp(1j*va); zr+1j*zi];
        end

        function names = opf_legacy_user_var_names(obj)
            names = {'Va', 'Vm', 'Pg', 'Qg'};
        end
    end     %% methods
end         %% classdef
