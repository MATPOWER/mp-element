classdef acsp_aggregate < ac_aggregate & acsp_model

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
        %% constructor
        function obj = acsp_aggregate(varargin)
            obj@ac_aggregate(varargin{:});
            obj.element_classes = ...
                { @acsp_bus, @ac_gen, @ac_load, @ac_branch };
            if isempty(obj.node)    %% skip if constructed from existing object
                obj.init_set_types();   %% should be called in mp_idx_manager
                                        %% constructor, if not for:
                                        %% https://savannah.gnu.org/bugs/?52614
            end
        end

        function obj = def_set_types(obj)
            def_set_types@ac_aggregate(obj);        %% call parent first
            obj.set_types.va = 'VOLTAGE ANG VARS (va)';
            obj.set_types.vm = 'VOLTAGE MAG VARS (vm)';
        end

        function d2G = opf_power_balance_hess(obj, x, lam)
            nlam = length(lam) / 2;
            lamP = lam(1:nlam);
            lamQ = lam((1:nlam)+nlam);

            d2Gr = obj.port_inj_power_hess(x, lamP);
            d2Gi = obj.port_inj_power_hess(x, lamQ);

            d2G = real(d2Gr) + imag(d2Gi);
        end
    end     %% methods
end         %% classdef
