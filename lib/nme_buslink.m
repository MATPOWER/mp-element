classdef nme_buslink < nm_element & mp_form_acp

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
            obj.N = [ repmat(speye(obj.nk), 1, obj.nz);
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

        function obj = pf_add_constraints(obj, mm, nm, dm, mpopt)
            %% form constraint matrices for matching
            %%  voltage angles for pv and pq buses
            %%  voltage magnitudes for pq buses
            
            ad = mm.aux_data;

            %% voltage angle constraints
            A1 = ( obj.C([ad.pv; ad.pq], :) * obj.N )';
            k1 = find(any(A1, 2));
            n1 = length(k1);
            ang120 = 2*pi/3*ones(obj.nk, 1);
            b1 = [zeros(obj.nk, 1); ang120; -ang120];
            b1 = b1(k1);
            vs1 = struct('name', {'Va_pv', 'Va3_pv', 'Va3_pv', 'Va3_pv', ...
                                    'Va_pq', 'Va3_pq', 'Va3_pq', 'Va3_pq'}, ...
                            'idx', {{}, {1}, {2}, {3}, {}, {1}, {2}, {3}} );
            mm.add_lin_constraint('buslink_va', A1(k1, :), b1, b1, vs1);

            %% voltage magnitude constraints
            A2 = ( obj.C(ad.pq, :) * obj.N )';
            k2 = find(any(A2, 2));
            n2 = length(k2);
            b2 = zeros(n2, 1);
            vs2 = struct('name', {'Vm_pq', 'Vm3_pq', 'Vm3_pq', 'Vm3_pq'}, ...
                            'idx', {{}, {1}, {2}, {3}} );
            mm.add_lin_constraint('buslink_vm', A2(k2, :), b2, b2, vs2);
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
