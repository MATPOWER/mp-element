classdef mme_bus_nld_opf_acps_node_test < mme_bus_opf_acp

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        %% constructor
        function obj = mme_bus_nld_opf_acps_node_test()
            obj@mme_bus_opf_acp();
            obj.name = 'bus_nld';
        end

        function x0 = opf_interior_x0(obj, mm, nm, dm, x0)
            varef1 = mm.opf_interior_va(nm, dm);
            vm = obj.opf_interior_vm(mm, nm, dm);
            vv = mm.get_idx();
            x0(vv.i1.(['va_' obj.name]):vv.iN.(['va_' obj.name])) = varef1; %% angles set to 1st ref angle
            x0(vv.i1.(['vm_' obj.name]):vv.iN.(['vm_' obj.name])) = vm;     %% voltage magnitudes
        end
    end     %% methods
end         %% classdef
