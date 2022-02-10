classdef mm_shared_opf_legacy < handle
%MM_SHARED_OPF_LEGACY  MATPOWER mathematical model for optimal power flow (OPF) problem.
%   ?
%
%   MM_SHARED_OPF_LEGACY ... optimal power flow ...
%
%   Properties
%       ? - ?
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        cost = [];
        mpc = struct();
    end     %% properties

    methods
        function obj = def_set_types_legacy(obj)
            obj.set_types = struct(...
                    'var',  'VARIABLES', ...
                    'lin',  'LINEAR CONSTRAINTS', ...
                    'nle',  'NONLIN EQ CONSTRAINTS', ...
                    'nli',  'NONLIN INEQ CONSTRAINTS', ...
                    'qdc',  'QUADRATIC COSTS', ...
                    'nlc',  'GEN NONLIN COSTS', ...
                    'cost', 'LEGACY COSTS'  ...
                );
        end

        function obj = init_set_types_legacy(obj)
            %% finish initializing data structures for each type
            es = struct();  %% empty struct
            obj.cost.data = struct( ...
                'N', es, ...
                'H', es, ...
                'Cw', es, ...
                'dd', es, ...
                'rh', es, ...
                'kk', es, ...
                'mm', es, ...
                'vs', es );
            obj.cost.params = [];
        end

        function obj = build_legacy(obj, nm, dm, mpopt)
            if strcmp(obj.form_tag, 'dc') && toggle_softlims(obj.mpc, 'status')
                %% user data required by toggle_softlims
                branch_nme = nm.elements.branch;
                [Bbr, pbr] = branch_nme.get_params(1:branch_nme.nk, {'B', 'p'});
                obj.userdata.Bf = Bbr * branch_nme.C';
                obj.userdata.Pfinj = pbr;
            end

            %% execute userfcn callbacks for 'formulation' stage
            if isfield(obj.mpc, 'userfcn')
                userfcn = obj.mpc.userfcn;
            else
                userfcn = [];
            end
            obj = run_userfcn(userfcn, 'formulation', obj, mpopt);
        end

        function add_legacy_user_vars(obj, nm, dm, mpopt)
            %% save data
            obj.userdata.user_vars = obj.legacy_user_var_names();

            %% add any user-defined vars
            if isfield(dm.userdata.legacy_opf_user_mods, 'z')
                z = dm.userdata.legacy_opf_user_mods.z;
                if z.nz > 0
                    obj.add_var('z', z.nz, z.z0, z.zl, z.zu);
                    obj.userdata.user_vars{end+1} = 'z';
                end
            end
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

        function add_legacy_user_constraints_ac(obj, nm, dm, mpopt)
            obj.add_legacy_user_constraints(nm, dm, mpopt);

            if ~isempty(dm.userdata.legacy_opf_user_mods)
                uc = dm.userdata.legacy_opf_user_mods.nlc;
                for k = 1:length(uc)
                    obj.add_nln_constraint(uc{k}{:});
                end
            end
        end

    end     %% methods
end         %% classdef
