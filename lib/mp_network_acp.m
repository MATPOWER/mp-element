classdef mp_network_acp < mp_network_ac% & mp_form_acp

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
                    @nme_branch_acp, @nme_shunt_acp, ...
                    @nme_bus3p_acp, @nme_gen3p_acp, ...
                    @nme_load3p, @nme_line3p, @nme_xfmr3p, @nme_buslink_acp };

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

        function va = initial_voltage_angle(obj, idx)
            va = obj.params_var('va');
            if nargin > 1 && ~isempty(idx)
                va = va(idx);
            end
        end
    end     %% methods
end         %% classdef
