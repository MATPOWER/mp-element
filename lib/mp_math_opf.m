classdef mp_math_opf < mp_math
%MP_MATH_OPF  MATPOWER mathematical model for optimal power flow (OPF) problem.
%   ?
%
%   MP_MATH_OPF ... optimal power flow ...
%
%   Properties
%       ? - ?
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = build(obj, nm, dm, mpopt)
            build@mp_math(obj, nm, dm, mpopt);  %% call parent
            obj.add_aux_data(nm, dm, mpopt);
            obj.add_vars(nm, dm, mpopt);
            obj.add_constraints(nm, dm, mpopt);
            obj.add_costs(nm, dm, mpopt);
        end

        function obj = add_vars(obj, nm, dm, mpopt)
            add_vars@mp_math(obj, nm, dm, mpopt);   %% call parent

            %% legacy user-defined variables
            if isfield(dm.userdata, 'legacy_opf_user_mods')
                obj.add_legacy_user_vars(nm, dm, mpopt);
            end
        end

        function obj = add_system_vars(obj, nm, dm, mpopt)
            %% add network voltage and non-voltage state variables
            vars = horzcat(nm.model_vvars(), nm.model_zvars());
            for vtype = vars
                st = nm.(vtype{1});     %% set type
                d = st.data;
                mmx_i1 = obj.var.N + 1;
                for k = 1:st.NS
                    name = st.order(k).name;
                    idx = st.order(k).idx;
                    if isempty(idx)
                        obj.add_var(name, st.idx.N.(name), d.v0.(name), d.vl.(name), d.vu.(name), d.vt.(name));
                    else
                        if all(cell2mat(idx) == 1)
                            dim = size(st.idx.N.(name));
                            if dim(end) == 1, dim(end) = []; end   %% delete trailing 1
                            obj.init_indexed_name('var', name, num2cell(dim));
                        end
                        sn = struct('type', {'()'}, 'subs', {idx});
                        sc = sn;
                        sc.type = '{}';
                        N = subsref(st.idx.N.(name), sn);
                        v0 = subsref(d.v0.(name), sc);
                        vl = subsref(d.vl.(name), sc);
                        vu = subsref(d.vu.(name), sc);
                        vt = subsref(d.vt.(name), sc);
                        obj.add_var(name, idx, N, v0, vl, vu, vt);
                    end
                end
                mmx_iN = obj.var.N;
                obj.aux_data.var_map{end+1} = ...
                    {vtype{1}, [], [], [], mmx_i1, mmx_iN, []};
            end
        end

        function add_legacy_user_vars(obj, nm, dm, mpopt)
            %% save data
            obj.userdata.user_vars = obj.opf_legacy_user_var_names();

            %% add any user-defined vars
            if isfield(dm.userdata.legacy_opf_user_mods, 'z')
                z = dm.userdata.legacy_opf_user_mods.z;
                if z.nz > 0
                    obj.add_var('z', z.nz, z.z0, z.zl, z.zu);
                    obj.userdata.user_vars{end+1} = 'z';
                end
            end
        end

        function obj = add_system_constraints(obj, nm, dm, mpopt)
            %% call parent
            add_system_constraints@mp_math(obj, nm, dm, mpopt);

            %% legacy user-defined constraints
            if isfield(dm.userdata, 'legacy_opf_user_mods')
                obj.add_legacy_user_constraints(nm, dm, mpopt);
            end
        end

        function add_legacy_user_constraints(obj, nm, dm, mpopt)
            %% user-defined linear constraints
            if isfield(dm.userdata.legacy_opf_user_mods, 'lin')
                lin = dm.userdata.legacy_opf_user_mods.lin;
                if lin.nlin
                    uv = obj.get_userdata('user_vars');
                    obj.add_lin_constraint('usr', lin.A, lin.l, lin.u, uv);
                end
            end
        end

        function obj = add_system_costs(obj, nm, dm, mpopt)
        end

        function add_legacy_user_costs(obj, nm, dm, dc)
            if isfield(dm.userdata.legacy_opf_user_mods, 'cost') && ...
                    dm.userdata.legacy_opf_user_mods.cost.nw
                user_cost = dm.userdata.legacy_opf_user_mods.cost;
                uv = obj.get_userdata('user_vars');
                obj.add_legacy_cost('usr', user_cost, uv);

                %% implement legacy user costs using quadratic or general non-linear costs
                cp = obj.params_legacy_cost();  %% construct/fetch the parameters
                [N, H, Cw, rh, m] = deal(cp.N, cp.H, cp.Cw, cp.rh, cp.mm);
                [nw, nx] = size(N);
                if nw
                    if any(cp.dd ~= 1) || any(cp.kk)    %% not simple quadratic form
                        if dc                           %% (includes "dead zone" or
                            if any(cp.dd ~= 1)          %%  quadratic "penalty")
                                error('mp_network/add_legacy_user_costs: DC OPF can only handle legacy user-defined costs with d = 1');
                            end
                            if any(cp.kk)
                                error('mp_network/add_legacy_user_costs: DC OPF can only handle legacy user-defined costs with no "dead zone", i.e. k = 0');
                            end
                        else
                            %% use general nonlinear cost to implement legacy user cost
                            legacy_cost_fcn = @(x)opf_legacy_user_cost_fcn(x, cp);
                            obj.add_nln_cost('usr', 1, legacy_cost_fcn);
                        end
                    else                                %% simple quadratic form
                        %% use a quadratic cost to implement legacy user cost
                        %% f = 1/2 * w'*H*w + Cw'*w, where w = diag(m)*(N*x - rh)
                        %% Let: MN = diag(m)*N
                        %%      MR = M * rh
                        %%      HMR  = H  * MR;
                        %%      HtMR = H' * MR;
                        %%  =>   w = MN*x - MR
                        %% f = 1/2 * (MN*x - MR)'*H*(MN*x - MR) + Cw'*(MN*x - MR)
                        %%   = 1/2 * x'*MN'*H*MN*x +
                        %%          (Cw'*MN - 1/2 * MR'*(H+H')*MN)*x +
                        %%          1/2 * MR'*H*MR - Cw'*MR
                        %%   = 1/2 * x'*Q*w + c'*x + k

                        M    = sparse(1:nw, 1:nw, m, nw, nw);
                        MN   = M * N;
                        MR   = M * rh;
                        HMR  = H  * MR;
                        HtMR = H' * MR;
                        Q = MN' * H * MN;
                        c = full(MN' * (Cw - 1/2*(HMR+HtMR)));
                        k = (1/2 * HtMR - Cw)' * MR;
                        obj.add_quad_cost('usr', Q, c, k);
                    end
                end
            end
        end

        function varef1 = opf_interior_va(obj, nm, dm)
            %% return scalar va equal to angle of first reference node
            ad = obj.aux_data;
            ref1 = ad.ref(1);
            varef1 = ad.va(ref1);
        end

        function dm = data_model_update(obj, nm, dm, mpopt)
            nm.opf_data_model_update(obj, nm, dm, mpopt);
        end

        function nm = network_model_x_soln(obj, nm)
            %% convert solved state from math model to network model soln
            [nm.soln.v, nm.soln.z, nm.soln.x] = ...
                obj.opf_convert_x(obj.soln.x, nm, obj.aux_data);
        end
    end     %% methods
end         %% classdef
