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

        function [A, b_va, b_vm, Istack_] = voltage_constraints(obj)
            %% form constraint matrices for matching voltages
            nk = obj.nk;

            %% basic constraint matrix for voltage equality (all nodes)
            Istack = repmat(speye(nk), obj.nz, 1);  %% stacked identities
            A = [ Istack -speye(nk * obj.nz) ] * obj.C';
            ang120 = 2*pi/3*ones(nk, 1);

            %% RHS
            b_va = [zeros(nk, 1); ang120; -ang120]; %% angle constraints
            b_vm = zeros(nk*obj.nz, 1);             %% magnitude constraints

            if nargout > 3
                Istack_ = Istack;
            end
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

        function obj = pf_add_constraints(obj, mm, nm, dm, mpopt)
            %% add Qlink constraints for PV node on 3-phase side

            %% indices of buslinks connected to PV nodes via port 2 (3&4)
            ad = mm.aux_data;
            [~, k, ~] = find(obj.C(ad.pv, obj.nk+1:2*obj.nk));
            n = length(k);
            if n
                I = sparse(1:n, k, 1, n, obj.nk);
                zz = sparse(n, obj.nk);
                A = [I -I zz; I zz -I];
                b = zeros(2*n, 1);
                vs = struct('name', {'Qlink', 'Qlink', 'Qlink'}, ...
                            'idx', {{1}, {2}, {3}} );
                mm.add_lin_constraint('buslink_qlink', A, b, b, vs);
            end
        end

        function [A_va_pq, A_va_pv, b_va, A_vm, b_vm] = pf_voltage_constraints(obj, ad)
            %% form constraint matrices for matching
            %%  voltage angles for pv and pq nodes
            %%  voltage magnitudes for pq nodes
            %% columns of A_va_pv correspond to ...
            %%  'Va_pv', 'Va3_pv', 'Va3_pv', 'Va3_pv',
            %% columns of A_va_pq correspond to ...
            %%  'Va_pq', 'Va3_pq', 'Va3_pq', 'Va3_pq'
            %% columns of A_vm correspond to ...
            %%  'Vm_pq', 'Vm3_pq', 'Vm3_pq', 'Vm3_pq'
            nk = obj.nk;

            %% basic constraint matrix for voltage equality (all nodes)
            [A, b_va, b_vm, Istack] = obj.voltage_constraints();

            %% sub-matrices by node type
            A_ref = A(:, ad.ref);
            A_pv = A(:, ad.pv);
            A_pq = A(:, ad.pq);

            %% voltage angle constraints (all buslinks)
            A_va_pq = A_pq;
            A_va_pv = A_pv;

            %% voltage magnitude constraints (all buslinks)
            A_vm = A_pq;

            %% indices of buslinks connected to REF nodes (fixed va)
            [~, k_va1, ~] = find(obj.C(ad.ref, 1:nk));      %% via port 1
            [~, k_va2, ~] = find(obj.C(ad.ref, nk+1:2*nk)); %% via port 2 (3&4)

            %% adjust RHS of va constraints involving fixed voltage angles
            if ~isempty(k_va1) || ~isempty(k_va2)
                [va, vm] = obj.aux_data_va_vm(ad);
                va1 = obj.C(:, 1:nk)' * va;         %% port 1 voltage angles
                va2 = obj.C(:, nk+1:2*nk)' * va;    %% port 2 voltage angles
                va_adj = zeros(nk, 1);
                va_adj(k_va1) = -va1(k_va1);
                va_adj(k_va2) =  va2(k_va2);
                b_va = b_va + Istack * va_adj;

                %% delete redundate (and currently incorrect) constraints
                %% corresponding to k_va2 and ports 3 and 4
                if ~isempty(k_va2)
                    d_va = [nk+k_va2; 2*nk+k_va2];
                    A_va_pq(d_va, :) = [];
                    A_va_pv(d_va, :) = [];
                    b_va(d_va) = [];
                end
            end

            %% indices of buslinks connected to REF/PV nodes (fixed vm)
            rpv = [ad.ref;ad.pv];   %% indices of REF/PV nodes
            [~, k_vm1, ~] = find(obj.C(rpv, 1:nk));         %% via port 1
            [~, k_vm2, ~] = find(obj.C(rpv, nk+1:2*nk));    %% via port 2 (3&4)

            %% adjust RHS of vm constraints involving fixed voltage magnitudes
            if ~isempty(k_vm1) || ~isempty(k_vm2)
                [va, vm] = obj.aux_data_va_vm(ad);
                vm1 = obj.C(:, 1:nk)' * vm;         %% port 1 voltage magnitudes
                vm2 = obj.C(:, nk+1:2*nk)' * vm;    %% port 2 voltage magnitudes
                vm_adj = zeros(nk, 1);
                vm_adj(k_vm1) = -vm1(k_vm1);
                vm_adj(k_vm2) =  vm2(k_vm2);
                b_vm = b_vm + Istack * vm_adj;

                %% delete redundate (and currently incorrect) constraints
                %% corresponding to k_vm2 and ports 3 and 4
                if ~isempty(k_vm2)
                    d_vm = [nk+k_vm2; 2*nk+k_vm2];
                    A_vm(d_vm, :) = [];
                    b_vm(d_vm) = [];
                end
            end
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

        %%-----  CPF methods  -----
        function obj = cpf_add_constraints(obj, mm, nm, dm, mpopt)
            obj.pf_add_constraints(mm, nm, dm, mpopt);
        end
    end     %% methods
end         %% classdef
