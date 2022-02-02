classdef mme_bus_opf_acp < mme_bus_opf_ac

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'bus';
%     end
    
    methods
        function x0 = opf_interior_x0(obj, mm, nm, dm, x0)
            vv = mm.get_idx();
            varef1 = mm.opf_interior_va(nm, dm);
            vm = obj.opf_interior_vm(mm, nm, dm);
            x0(vv.i1.Va:vv.iN.Va) = varef1; %% angles set to 1st ref angle
            x0(vv.i1.Vm:vv.iN.Vm) = vm;     %% voltage magnitudes
        end
    end     %% methods
end         %% classdef
