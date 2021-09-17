classdef nme_buslink < nm_element %& mp_form_ac

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        %% constructor
        function obj = nme_buslink()
            obj@nm_element();
            obj.name = 'buslink';
            obj.np = 4;             %% this is a 4 port element
            obj.nz = 3;
        end

        function obj = add_zvars(obj, nm, dm, idx)
            p = idx{1};
            ng = obj.nk;

            if p == 1
                nm.init_indexed_name('zr', 'Plink', {obj.nz});
                nm.init_indexed_name('zi', 'Qlink', {obj.nz});
            end
            nm.add_var('zr', 'Plink', {p}, obj.nk, 0, -Inf, Inf);
            nm.add_var('zi', 'Qlink', {p}, obj.nk, 0, -Inf, Inf);
        end

        function obj = build_params(obj, nm, dm)
            build_params@nm_element(obj, nm, dm);   %% call parent
            I = (dm.base_kva / dm.base_mva / 1000) * speye(obj.nk);
            obj.N = [ repmat(I, 1, obj.nz);
                     -speye(obj.nk * obj.nz) ];
        end

        %%-----  PF methods  -----
        function obj = pf_add_vars(obj, mm, nm, dm, mpopt)
            mm.init_indexed_name('var', 'Plink', {3});
            mm.init_indexed_name('var', 'Qlink', {3});
            mmx_i1 = mm.var.N + 1;
            for p = 1:3
                mm.add_var('Plink', {p}, obj.nk, 0, -Inf, Inf);
            end
            mmx_iN = mm.var.N;
            mm.aux_data.var_map{end+1} = ...
                {'zr', nm.zr.idx.i1.Plink(1), nm.zr.idx.iN.Plink(end), [], mmx_i1, mmx_iN, []};

            mmx_i1 = mm.var.N + 1;
            for p = 1:3
                mm.add_var('Qlink', {p}, obj.nk, 0, -Inf, Inf);
            end
            mmx_iN = mm.var.N;
            mm.aux_data.var_map{end+1} = ...
                {'zi', nm.zi.idx.i1.Qlink(1), nm.zi.idx.iN.Qlink(end), [], mmx_i1, mmx_iN, []};
        end

        function [A_va, b_va, A_vm, b_vm] = pf_voltage_constraints(obj, ad)
            %% form constraint matrices for matching
            %%  voltage angles for pv and pq nodes
            %%  voltage magnitudes for pq nodes
            %% columns of A_va correspond to ...
            %%  'Va_pv', 'Va3_pv', 'Va3_pv', 'Va3_pv',
            %%  'Va_pq', 'Va3_pq', 'Va3_pq', 'Va3_pq'
            %% columns of A_vm correspond to ...
            %%  'Vm_pq', 'Vm3_pq', 'Vm3_pq', 'Vm3_pq'

            A = [ repmat(speye(obj.nk), obj.nz, 1) ...
                    -speye(obj.nk * obj.nz) ] * obj.C';
            ang120 = 2*pi/3*ones(obj.nk, 1);

            %% voltage angle constraints
            A_va = A(:, [ad.pv; ad.pq]);
            k_va = find(any(A_va, 2));
            b_va = [zeros(obj.nk, 1); ang120; -ang120];
            A_va = A_va(k_va, :);
            b_va = b_va(k_va);

            %% voltage magnitude constraints
            A_vm = A(:, ad.pq);
            k_vm = find(any(A_vm, 2));
            b_vm = zeros(length(k_vm), 1);
            A_vm = A_vm(k_vm, :);
        end

%         function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
%             dme = obj.data_model_element(dm);
%             ss = nm.get_idx('state');
% 
%             for p = 1:obj.nz
%                 %% generator active power
%                 sg = nm.soln.z(ss.i1.buslink(p):ss.iN.buslink(p)) * dm.base_kva;
%                 pg = real(sg);
% 
%                 %% update in the data model
%                 dme.tab.(sprintf('pg%d', p))(dme.on) = real(sg);
%                 dme.tab.(sprintf('pf%d', p))(dme.on) = cos(angle(sg));
%             end
%         end
    end     %% methods
end         %% classdef
