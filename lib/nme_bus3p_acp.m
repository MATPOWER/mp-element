classdef nme_bus3p_acp < nme_bus3p & mp_form_acp

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'bus3p';
%     end

    methods
        function obj = add_vvars(obj, nm, dm, idx)
            dme = obj.data_model_element(dm);
            nb = obj.nk;
            p = idx{1};

            %% prepare angle bounds for ref buses
            ref = dme.type == NODE_TYPE.REF;
            va_lb = -Inf(nb, 1);
            va_ub =  Inf(nb, 1);
            vm_start = dme.(sprintf('vm%d_start', p));
            va_start = dme.(sprintf('va%d_start', p));
            va1_lb(ref) = va_start(ref);
            va1_ub(ref) = va_start(ref);

            if p == 1
                nm.init_indexed_name('va', 'Va3', {obj.nn});
                nm.init_indexed_name('vm', 'Vm3', {obj.nn});
            end
            nm.add_var('va', 'Va3', {p}, nb, va_start, va_lb, va_ub);
            nm.add_var('vm', 'Vm3', {p}, nb, vm_start, 0, 2);
        end

        %%-----  OPF methods  -----
        function x0 = opf_interior_x0(obj, mm, nm, dm, x0)
            vv = mm.get_idx();
            vm = 1;     %% voltage magnitude set to 1 p.u.
            %% voltage angles set to angle of 1st ref node
            %% + phase offset, i.e. +/-120 deg
            varef1 = nm.opf_interior_va(mm, dm);

            x0(vv.i1.Va3(1):vv.iN.Va3(1)) = varef1;
            x0(vv.i1.Va3(2):vv.iN.Va3(2)) = varef1-2*pi/3;
            x0(vv.i1.Va3(3):vv.iN.Va3(3)) = varef1+2*pi/3;
            x0(vv.i1.Vm3(1):vv.iN.Vm3(1)) = vm;
            x0(vv.i1.Vm3(2):vv.iN.Vm3(2)) = vm;
            x0(vv.i1.Vm3(3):vv.iN.Vm3(3)) = vm;
        end
    end     %% methods
end         %% classdef
