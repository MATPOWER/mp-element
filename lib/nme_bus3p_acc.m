classdef nme_bus3p_acc < nme_bus3p & mp_form_acc

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

            %% prepare box bounds for voltage coordinates
            vm_start = dme.(sprintf('vm%d_start', p));
            va_start = dme.(sprintf('va%d_start', p));
            v_start = vm_start .* exp(1j * va_start);
            vclim = 1.5;

            if p == 1
                nm.init_indexed_name('vr', 'Vr3', {obj.nn});
                nm.init_indexed_name('vi', 'Vi3', {obj.nn});
            end
            nm.add_var('vr', 'Vr3', {p}, nb, real(v_start), -vclim, vclim);
            nm.add_var('vi', 'Vi3', {p}, nb, imag(v_start), -vclim, vclim);
        end

        %%-----  OPF methods  -----
        function x0 = opf_interior_x0(obj, mm, nm, dm, x0)
            vv = mm.get_idx();
            vm = 1;     %% voltage magnitude set to 1 p.u.
            %% voltage angles set to angle of 1st ref node
            %% + phase offset, i.e. +/-120 deg
            varef1 = nm.opf_interior_va(mm, dm);
            v1 = vm * exp(1j*varef1);
            v2 = vm * exp(1j*(varef1-2*pi/3));
            v3 = vm * exp(1j*(varef1+2*pi/3));

            x0(vv.i1.Vr3(1):vv.iN.Vr3(1)) = real(v1);
            x0(vv.i1.Vr3(2):vv.iN.Vr3(2)) = real(v2);
            x0(vv.i1.Vr3(3):vv.iN.Vr3(3)) = real(v3);
            x0(vv.i1.Vi3(1):vv.iN.Vi3(1)) = imag(v1);
            x0(vv.i1.Vi3(2):vv.iN.Vi3(2)) = imag(v2);
            x0(vv.i1.Vi3(3):vv.iN.Vi3(3)) = imag(v3);
        end
    end     %% methods
end         %% classdef
