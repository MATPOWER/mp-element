classdef mp_data_mpc2 < mp_data
%MP_DATA_MPC2  Implementation of MATPOWER data model for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpc
    end     %% properties

    methods
        %% constructor
        function obj = mp_data_mpc2(m, class_list)
            %% call parent constructor
            obj@mp_data();
            obj.element_classes = ...
                { @dme_bus_mpc2, @dme_gen_mpc2, @dme_load_mpc2, ...
                    @dme_branch_mpc2, @dme_shunt_mpc2 };

            if nargin
                if nargin > 1
                    obj.modify_element_classes(class_list);
                end
                %% load case and create mappings
                obj.mpc = loadcase(m);
                obj.build();
            end
        end

        function ref = node_type_ref(obj, node_type)
            %% define constants
            [PQ, PV, REF, NONE] = idx_bus;

            ref = find(node_type == REF);
        end

        function pv = node_type_pv(obj, node_type)
            %% define constants
            [PQ, PV, REF, NONE] = idx_bus;

            pv = find(node_type == PV);
        end

        function pq = node_type_pq(obj, node_type)
            %% define constants
            [PQ, PV, REF, NONE] = idx_bus;

            pq = find(node_type == PQ);
        end

        function obj = ext2int(obj, mpopt)
            if ~isfield(obj.mpc, 'order') || obj.mpc.order.state == 'e'
                if nargin > 1
                    obj.mpc = ext2int(obj.mpc, mpopt);
                else
                    obj.mpc = ext2int(obj.mpc);
                end
            else
%                 warning('mp_data_mpc2/ext2int: data model already in internal format');
            end
        end

        function obj = int2ext(obj, mpopt)
            if isfield(obj.mpc, 'order') && obj.mpc.order.state == 'i'
                if nargin > 1
                    obj.mpc = int2ext(obj.mpc, mpopt);
                else
                    obj.mpc = int2ext(obj.mpc);
                end
            else
%                 warning('mp_data_mpc2/int2ext: data model already in external format');
            end
        end

        function display(obj)
            fprintf('Data Model class : %s\n', class(obj));
            mpc = obj.mpc
        end

        function print_soln(obj, fname)
        end

        function save_soln(obj, fname)
        end

        %%-----  PF methods  -----
        function [dm1, dm2] = fdpf_B_matrix_models(obj, alg)
            %% [dmp, dmpp] = obj.fdpf_B_matrix_models(alg)
            %% dmpp = obj.fdpf_B_matrix_models(alg)
            %% returns copies of dm used for building B prime, B double prime
            %% for fast-decoupled power flow

            %% define named indices into bus, branch matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

            %% modify data model to form Bp (B prime)
            if nargout > 1      %% for both Bp and Bpp
                mpc1 = obj.mpc;
                mpc1.bus(:, BS) = 0;            %% zero out shunts at buses
                mpc2 = mpc1;
                mpc1.branch(:, BR_B) = 0;       %% zero out line charging shunts
                mpc1.branch(:, TAP) = 1;        %% cancel out taps
                if strcmp(alg, 'FDXB')          %% if XB method
                    mpc1.branch(:, BR_R) = 0;   %% zero out line resistance
                end
                dm1 = feval(class(obj), mpc1);
            else
                mpc2 = obj.mpc;
            end

            %% modify data model to form Bpp (B double prime)
            mpc2.branch(:, SHIFT) = 0;      %% zero out phase shifters
            if strcmp(alg, 'FDBX')          %% if BX method
                mpc2.branch(:, BR_R) = 0;   %% zero out line resistance
            end

            if nargout > 1      %% for both Bp and Bpp
                dm2 = feval(class(obj), mpc2);
            else                %% for just Bpp
                dm1 = feval(class(obj), mpc2);
            end
        end

        %%-----  OPF methods  -----
        function [A, l, u, i] = branch_angle_diff_constraint(obj, ignore);
            baseMVA = obj.mpc.baseMVA;
            branch = obj.elm_by_name('branch').get_table(obj);
            nb = obj.elm_by_name('bus').n;
            mpopt = struct('opf', struct('ignore_angle_lim', ignore));
            [A, l, u, i]  = makeAang(baseMVA, branch, nb, mpopt);
        end

        function [Ah, uh, Al, ul, data] = gen_pq_capability_constraint(obj);
            baseMVA = obj.mpc.baseMVA;
            gen_dme = obj.elm_by_name('gen');
            gen = gen_dme.get_table(obj);

            [Ah, uh, Al, ul, data] = makeApq(baseMVA, gen(gen_dme.on, :));
        end

        function [uv, z] = opf_legacy_user_vars(obj, uv_names, nx, mpopt)
            %% create (read-only) copies of individual fields for convenience
            mpc = obj.mpc;
            [baseMVA, bus, gen, branch, gencost, Au, lbu, ubu, mpopt, ...
                N, fparm, H, Cw, z0, zl, zu, userfcn] = opf_args(mpc, mpopt);

            %% get some more problem dimensions
            if isfield(mpc, 'A')
                nlin = size(mpc.A, 1);  %% number of linear user constraints
            else
                nlin = 0;
            end
            if isfield(mpc, 'N')
                nw = size(mpc.N, 1);    %% number of general cost vars, w
            else
                nw = 0;
            end

            %% get number of user vars, check consistency
            if nlin
                nz = size(mpc.A, 2) - nx; %% number of user z variables
                if nz < 0
                    error('mp_data_mpc2/add_opf_legacy_user_vars: user supplied A matrix must have at least %d columns.', nx);
                end
            else
                nz = 0;               %% number of user z variables
                if nw                 %% still need to check number of columns of N
                    if size(mpc.N, 2) ~= nx;
                        error('mp_data_mpc2/add_opf_legacy_user_vars: user supplied N matrix must have %d columns.', nx);
                    end
                end
            end

            %% package up return data
            uv = struct( ...
                    'user_vars', {uv_names}, ...
                    'nlin', nlin, ...
                    'nw', nw, ...
                    'A', Au, ...
                    'l', lbu, ...
                    'u', ubu, ...
                    'N', N, ...
                    'fparm', fparm, ...
                    'H', H, ...
                    'Cw', Cw ...
                );
            z = struct('nz', nz, 'z0', z0, 'zl', zl, 'zu', zu);
        end

        function uc = opf_legacy_user_constraints(obj, uv_names, nx, mpopt)
            mpc = dm.mpc;

            %% check for user-defined nonlinear constraints
            nnle = 0;   %% number of nonlinear user-defined equality cons
            nnli = 0;   %% number of nonlinear user-defined inequality cons
            if isfield(mpc, 'user_constraints')
                if isfield(mpc.user_constraints, 'nle')
                    for k = 1:length(mpc.user_constraints.nle)
                        nnle = nnle + mpc.user_constraints.nle{k}{2};
                    end
                end
                if isfield(mpc.user_constraints, 'nli')
                    for k = 1:length(mpc.user_constraints.nli)
                        nnli = nnli + mpc.user_constraints.nli{k}{2};
                    end
                end
            end

            %% initialize cell array for add_nln_constraint() args
            uc = cell(nnle+nnli, 1);

            %% user-defined nonlinear equalities
            if nnle
                for k = 1:length(mpc.user_constraints.nle)
                    nlc = mpc.user_constraints.nle{k};
                    fcn  = eval(['@(x)' nlc{3} '(x, nlc{6}{:})']);
                    hess = eval(['@(x, lam)' nlc{4} '(x, lam, nlc{6}{:})']);
                    uc{k} = {nlc{1:2}, 1, fcn, hess, nlc{5}};
                end
            end

            %% user-defined nonlinear inequalities
            if nnli
                for k = 1:length(mpc.user_constraints.nli)
                    nlc = mpc.user_constraints.nli{k};
                    fcn  = eval(['@(x)' nlc{3} '(x, nlc{6}{:})'])
                    hess = eval(['@(x, lam)' nlc{4} '(x, lam, nlc{6}{:})'])
                    uc{nnle+k} = {nlc{1:2}, 0, fcn, hess, nlc{5}};
                end
            end
        end

        function mm = run_userfcn(obj, mm, mpopt)
            %% execute userfcn callbacks for 'formulation' stage
            mpc = obj.mpc;
            if isfield(mpc, 'userfcn')
                userfcn = mpc.userfcn;
            else
                userfcn = [];
            end
            mm = run_userfcn(userfcn, 'formulation', mm, mpopt);
        end
    end     %% methods
end         %% classdef
