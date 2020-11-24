classdef mpe_network_acp < mpe_network_ac & mp_model_acp

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
        function obj = mpe_network_acp()
            obj@mpe_network_ac();
            obj.element_classes = ...
                { @mpe_bus_acp, @mpe_gen_acp, @mpe_load_acp, @mpe_branch_acp, @mpe_shunt_acp };

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  INIT_SET_TYPES() is called explicitly as needed
            %%              (if obj.node is empty) in CREATE_MODEL() and
            %%              DISPLAY(), after object construction, but before
            %%              object use.
        end

        function obj = def_set_types(obj)
            def_set_types@mpe_network_ac(obj);      %% call parent first
            obj.set_types.va = 'VOLTAGE ANG VARS (va)';
            obj.set_types.vm = 'VOLTAGE MAG VARS (vm)';
        end


        %%-----  OPF methods  -----
        function x_ = x2x_(obj, x)
            %% convert (real) opt_model x to (complex) network model x_
            nv_ = obj.nv / 2;       %% number of voltage vars (sysx=1)
            nz_ = obj.nz;           %% number of state vars
            va = x(1:nv_, :);       b = nv_;
            vm = x(b+1:b+nv_, :);   b = b + nv_;
            zr = x(b+1:b+nz_, :);   b = b + nz_;
            zi = x(b+1:b+nz_, :);
            x_ = [vm .* exp(1j*va); zr+1j*zi];
        end

        function names = opf_legacy_user_var_names(obj)
            names = {'Va', 'Vm', 'Pg', 'Qg'};
        end
    end     %% methods
end         %% classdef
