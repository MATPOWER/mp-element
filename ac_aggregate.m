classdef ac_aggregate < mp_aggregate% & mp_model_ac
%AC_AGGREGATE Abstract class, explicitly a subclass of mp_aggregate and
%             implicitly assumed to be subclasses of MP_MODEL_AC as well

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        zr = [];
        zi = [];
        inln_list = {};         %% private: list of indexes of mpe's w/inln
        snln_list = {};         %% private: list of indexes of mpe's w/snln
        inln_hess_list = {};    %% private: list of indexes of mpe's w/inln_hess
        snln_hess_list = {};    %% private: list of indexes of mpe's w/snln_hess
    end
    
    methods
        function obj = def_set_types(obj)
            def_set_types@mp_aggregate(obj);        %% call parent first
            obj.set_types.zr = 'NON-VOLTAGE VARS REAL (zr)';
            obj.set_types.zi = 'NON-VOLTAGE VARS IMAG (zi)';
        end

        function obj = build_params(obj, nm, mpc)
            %% call parent to build individual element parameters
            build_params@mp_aggregate(obj, nm, mpc);

            %% aggregate parameters from individual elements
            obj.Y = obj.stack_matrix_params('Y', 1);
            obj.L = obj.stack_matrix_params('L', 0);
            obj.M = obj.stack_matrix_params('M', 1);
            obj.N = obj.stack_matrix_params('N', 0);
            obj.i = obj.stack_vector_params('i');
            obj.s = obj.stack_vector_params('s');

            %% add general nonlinear function if any element has one defined
            for k = 1:length(obj.mpe_list)
                if ~isempty(obj.mpe_list{k}.inln)
                    obj.inln_list{end+1} = k;
                    if ~isempty(obj.mpe_list{k}.inln_hess)
                        obj.inln_hess_list{end+1} = k;
                    end
                end
                if ~isempty(obj.mpe_list{k}.snln)
                    obj.snln_list{end+1} = k;
                    if ~isempty(obj.mpe_list{k}.snln_hess)
                        obj.snln_hess_list{end+1} = k;
                    end
                end
            end
            if ~isempty(obj.inln_list)
                obj.inln = @(x_, sysx, idx)port_inj_nln(obj, 'i', x_, sysx, idx);
                if ~isempty(obj.inln_hess_list)
                    obj.inln_hess = @(x_, lam, sysx, idx)port_inj_nln_hess(obj, 'i', x_, lam, sysx, idx);
                end
            end
            if ~isempty(obj.snln_list)
                obj.snln = @(x_, sysx, idx)port_inj_nln(obj, 's', x_, sysx, idx);
                if ~isempty(obj.snln_hess_list)
                    obj.snln_hess = @(x_, lam, sysx, idx)port_inj_nln_hess(obj, 's', x_, lam, sysx, idx);
                end
            end
        end

        function [g, gv1, gv2, gzr, gzi] = port_inj_nln(obj, si, x_, sysx, idx)
            if nargin < 5
                idx = [];
                if nargin < 4
                    sysx = 1;
                end
            end

            %% current or power
            fcn = [si 'nln'];
            fcn_list = [fcn '_list'];

            %% initialize
            if isempty(idx)
                sel = 0;        %% all ports
                np = obj.np;
            else
                sel = 1;        %% selected ports only
                np = length(idx);
            end
            nv = obj.get_nv_(sysx);
            nz = obj.nz;
            nc = size(x_, 2);   %% num of cols in x_, for evaluating multiple x_
            g = zeros(np, nc);
            gv1 = sparse(np, nv);
            gv2 = sparse(np, nv);
            gzr = sparse(np, nz);
            gzi = sparse(np, nz);

            %% loop through elements w/gen nonlin fcns, evaluate them
            for kk = obj.(fcn_list)
                k = kk{1};      %% index into obj.mpe_list
                mpe = obj.mpe_list{k};
                i1 = obj.mpe_port_map(k, 1);    %% starting aggregate port index
                iN = obj.mpe_port_map(k, 2);    %% ending aggregate port index

                %% set up port index vector for mpe
                if sel
                    apidx = find(idx >= i1 & idx <= iN);   %% aggregate port indices in range
                    if isempty(apidx)  %% skip if selected ports, but none in range
                        continue;
                    end
                    mpe_idx = idx(apidx) - i1 + 1;  %% port index vector for mpe
                else
                    mpe_idx = [];                   %% all ports for mpe
                end

                %% set up proper x_ for mpe
                if sysx
                    mpe_x_ = x_;
                else
                    j1 = obj.mpe_z_map(k, 1);   %% starting aggregate z-var index
                    jN = obj.mpe_z_map(k, 2);   %% ending aggregate z-var index
                    mpe_x_ = [  x_(i1:iN, :);
                                x_(nv+j1:nv+jN, :)  ];
                end

                %% call nonlinear function
                gg = cell(1, nargout);
                [gg{:}] = mpe.(fcn)(mpe_x_, sysx, mpe_idx);

                %% insert the results in aggregate output args
                if sel
                    g(apidx, :) = gg{1};
                    if nargout > 1
                        if sysx
                            gv1(apidx, :) = gg{2};
                            gv2(apidx, :) = gg{3};
                            if nargout > 3 && mpe.nz
                                gzr(apidx, :) = gg{4};
                                gzi(apidx, :) = gg{5};
                            end
                        else
                            gv1(apidx, i1:iN) = gg{2};
                            gv2(apidx, i1:iN) = gg{3};
                            if nargout > 3 && mpe.nz
                                gzr(apidx, j1:jN) = gg{4};
                                gzi(apidx, j1:jN) = gg{5};
                            end
                        end
                    end
                else
                    g(i1:iN, :) = gg{1};
                    if nargout > 1
                        if sysx
                            gv1(i1:iN, :) = gg{2};
                            gv2(i1:iN, :) = gg{3};
                            if nargout > 3 && mpe.nz
                                gzr(i1:iN, :) = gg{4};
                                gzi(i1:iN, :) = gg{5};
                            end
                        else
                            gv1(i1:iN, i1:iN) = gg{2};
                            gv2(i1:iN, i1:iN) = gg{3};
                            if nargout > 3 && mpe.nz
                                gzr(i1:iN, j1:jN) = gg{4};
                                gzi(i1:iN, j1:jN) = gg{5};
                            end
                        end
                    end
                end
            end     %% for loop
        end

        function H = port_inj_nln_hess(obj, si, x_, lam, sysx, idx)
            % H = obj.port_inj_power_hess(x_, lam)
            % H = obj.port_inj_power_hess(x_, lam, sysx)
            % H = obj.port_inj_power_hess(x_, lam, sysx, idx)
            if nargin < 6
                idx = [];
                if nargin < 5
                    sysx = 1;
                end
            end

            %% current or power
            fcn = [si 'nln_hess'];
            fcn_list = [fcn '_list'];

            %% initialize
            n = 2 * length(x_);
            H = sparse(n, n);

            %% loop through elements w/gen nonlin Hessians, evaluate them
            for kk = obj.(fcn_list)
                k = kk{1};      %% index into obj.mpe_list
                mpe = obj.mpe_list{k};
                i1 = obj.mpe_port_map(k, 1);    %% starting aggregate port index
                iN = obj.mpe_port_map(k, 2);    %% ending aggregate port index

                %% set up x_ for mpe & corresp row/col indices for mpe
                if sysx
                    mpe_x_ = x_;
                else
                    nv = obj.get_nv_(sysx);
                    nz = obj.nz;
                    j1 = obj.mpe_z_map(k, 1);   %% starting aggregate z-var index
                    jN = obj.mpe_z_map(k, 2);   %% ending aggregate z-var index
                    mpe_x_ = [  x_(i1:iN, :);
                                x_(nv+j1:nv+jN, :)  ];

                    %% indices of rows/cols of H corresponding to mpe x_
                    h = [(i1:iN) nv+(i1:iN) 2*nv+(j1:jN) 2*nv+nz+(j1:jN)].';
                end
                
                %% set up port index and lambda vectors for mpe
                if ~isempty(idx)    %% selected ports only
                    apidx = find(idx >= i1 & idx <= iN);    %% aggregate port indices in range
                    if isempty(apidx)  %% skip if selected ports, but none in range
                        continue;
                    end
                    mpe_idx = idx(apidx) - i1 + 1;  %% port index vector for mpe
                    mpe_lam = lam(apidx);           %% corresponding lam
                else                %% all ports
                    mpe_idx = [];
                    mpe_lam = lam(i1:iN);
                end

                %% call nonlinear function
                mpe_H = mpe.(fcn)(mpe_x_, mpe_lam, sysx, mpe_idx);

                %% accumulate output
                if sysx
                    H = H + mpe_H;
                else
                    H(h,h) = H(h,h) + mpe_H;
                end
            end
        end


        %%-----  PF methods  -----
        function ntv = power_flow_node_types(obj, nm, mpc, idx)
%         function [ntv, nts] = power_flow_node_types(obj, nm, mpc, idx)
            %% create empty cell array for node type vectors
            tt = cell(length(obj.node.order), 1);
            
            %% get node type vector from each node-creating MPE
            for k = 1:length(obj.node.order)
                mpe = obj.mpe_by_name(obj.node.order(k).name);
                tt{k} = mpe.power_flow_node_types(obj, mpc, obj.state.order(k).idx);
            end

            %% concatenate into a single node type vector
            ntv = vertcat(tt{:});

%             %% create node type struct
%             if nargout > 1
%                 %% define constants
%                 [PQ, PV, REF, NONE] = idx_bus;
% 
%                 ref = find(ntv == REF);     %% reference node indices
%                 pv  = find(ntv == PV );     %% PV node indices
%                 pq  = find(ntv == PQ );     %% PQ node indices
%                 nts = struct('ref', ref, 'pv', pv, 'pq', pq, ...
%                     'nref', length(ref), 'npv', length(pv), 'npq', length(pq));
%             end
        end

        function ad = power_flow_aux_data(obj, mpc, mpopt)
            %% get model variables
            vvars = obj.model_vvars();
            zvars = obj.model_zvars();
            v1 = obj.params_var(vvars{1});
            v2 = obj.params_var(vvars{2});
            zr = obj.params_var(zvars{1});
            zi = obj.params_var(zvars{2});

            %% define constants
            [PQ, PV, REF, NONE] = idx_bus;

            %% get node types
            ntv = obj.power_flow_node_types(obj, mpc);
            ref = find(ntv == REF);     %% reference node indices
            pv  = find(ntv == PV );     %% PV node indices
            pq  = find(ntv == PQ );     %% PQ node indices

            %% create aux_data struct
            ad = struct( ...
                'v1', v1, ...               %% initial value of v1 (va or vr)
                'v2', v2, ...               %% initial value of v2 (vm or vi)
                'zr', zr, ...               %% initial value of zr
                'zi', zi, ...               %% initial value of zi
                'ref',  ref, ...            %% REF node indices
                'nref', length(ref), ...    %% number of REF nodes
                'pv',  pv, ...              %% PV node indices
                'npv', length(pv), ...      %% number of PV nodes
                'pq',  pq, ...              %% PQ node indices
                'npq', length(pq) ...       %% number of PQ nodes
            );
        end

        function add_pf_constraints(obj, nm, om, mpc, mpopt)
            %% system constraints
            obj.add_pf_system_constraints(om, mpc, mpopt);
            
%             %% each element adds its PF constraints
%             for mpe = obj.mpe_list
%                 mpe{1}.add_pf_constraints(nm, om, mpc, mpopt);
%             end
        end

        function add_pf_system_constraints(obj, om, mpc, mpopt)
            %% can be overridden to add additional system constraints

            %% node balance constraints
            obj.add_pf_node_balance_constraints(om, mpc, mpopt);
        end

        function om = setup_power_flow(obj, mpc, mpopt)
            %% MATPOWER options
            if nargin < 3
                mpopt = mpoption;
            end

            %% create mathematical model
            om = opt_model();

            %% construct power flow auxiliary data, save in MP-Opt-Model
            om.userdata.power_flow_aux_data = ...
                obj.power_flow_aux_data(mpc, mpopt);

            %% add variables and constraints
            if obj.np ~= 0      %% skip for empty model
                obj.add_pf_vars(obj, om, mpc, mpopt);
                obj.add_pf_constraints(obj, om, mpc, mpopt);
            end
        end

        function [v_, success, i, data] = solve_power_flow(obj, mpc, mpopt)
            %% MATPOWER options
            if nargin < 3
                mpopt = mpoption;
            end

            if mpopt.verbose, fprintf('-----  solve_power_flow()  -----\n'); end

            %% create MP-Opt-Model object
            om = obj.setup_power_flow(mpc, mpopt);

            %% solve it
            alg = 'DEFAULT';
            switch mpopt.pf.alg
                case 'DEFAULT'
                    opt = mpopt2nleqopt(mpopt, om.problem_type(), 'DEFAULT');
                case {'NR', 'NR-SP', 'NR-SC', 'NR-SH', 'NR-IP', 'NR-IC', 'NR-IH'}
                    opt = mpopt2nleqopt(mpopt, om.problem_type(), 'NEWTON');
                case {'FDXB', 'FDBX'}
                    opt = mpopt2nleqopt(mpopt, om.problem_type(), 'FD');
                    opt.fd_opt.jac_approx_fcn = @()obj.fd_jac_approx(om, mpc, mpopt);
                    opt.fd_opt.labels = {'P', 'Q'};
                case 'FSOLVE'
                    opt = mpopt2nleqopt(mpopt, om.problem_type(), 'FSOLVE');
                case 'GS'
                    opt = mpopt2nleqopt(mpopt, om.problem_type(), 'GS');
                    opt.gs_opt.x_update_fcn = ...
                        @(x, f)obj.gs_x_update(x, f, om, mpc, mpopt);
                otherwise
                    error('ac_aggregate/solve_power_flow: invalid value for MPOPT.PF.ALG (%s)', mpopt.pf.alg);
            end
            opt.verbose = mpopt.verbose;

%             [x, f, eflag, output, J] = om.solve(opt);
            [x, f, eflag, output] = om.solve(opt);
            success = (eflag > 0);

% om
% vv = om.get_idx('var');
% x(vv.i1.Pg:vv.iN.Pg) * mpc.baseMVA
% keyboard

            if isfield(output, 'iterations')
                i = output.iterations;
            else
                i = -1;
            end

            %% convert back to complex voltage vector
            [v_, z_] = obj.pfx2vz(x, om.get_userdata('power_flow_aux_data'));
        end


        %%-----  OPF methods  -----
        function [G, Gv1, Gv2, Gzr, Gzi] = nodal_complex_current_balance(obj, x_)
            %% node incidence matrix
            C = obj.C;

            %% get port current injections with derivatives
            if nargout > 1
                [I, Iv1, Iv2, Izr, Izi] = obj.port_inj_current(x_, 1);
                Gv1 = C * Iv1;      %% Gva or Gvr
                Gv2 = C * Iv2;      %% Gvm or Gvi
                Gzr = C * Izr;
                Gzi = C * Izi;
            else
                I = obj.port_inj_current(x_, 1);
            end

            %% nodal current balance
            G = C * I;
        end

        function [G, Gv1, Gv2, Gzr, Gzi] = nodal_complex_power_balance(obj, x_)
            %% node incidence matrix
            C = obj.C;

            %% get port power injections with derivatives
            if nargout > 1
                [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1);
                Gv1 = C * Sv1;      %% Gva or Gvr
                Gv2 = C * Sv2;      %% Gvm or Gvi
                Gzr = C * Szr;
                Gzi = C * Szi;
            else
                S = obj.port_inj_power(x_, 1);
            end

            %% nodal power balance
            G = C * S;
        end

        function d2G = nodal_complex_current_balance_hess(obj, x_, lam)
            %% get port power injection hessians
            d2G = obj.port_inj_current_hess(x_, obj.C' * lam);
        end

        function d2G = nodal_complex_power_balance_hess(obj, x_, lam)
            %% get port power injection hessians
            d2G = obj.port_inj_power_hess(x_, obj.C' * lam);
        end

        function [g, dg] = opf_current_balance_fcn(obj, x)
            x_ = obj.x2x_(x);           %% convert real to complex x
            if nargout > 1
                [G, Gv1, Gv2, Gzr, Gzi] = obj.nodal_complex_current_balance(x_);
                Gx = [Gv1 Gv2 Gzr Gzi];
                dg = [  real(Gx);       %% Re{I} mismatch w.r.t v1, v2, zr, zi
                        imag(Gx)    ];  %% Im{I} mismatch w.r.t v1, v2, zr, zi
            else
                G = obj.nodal_complex_current_balance(x_);
            end
            g = [ real(G);              %% real current mismatch
                  imag(G) ];            %% imaginary current mismatch
        end

        function [g, dg] = opf_power_balance_fcn(obj, x)
            x_ = obj.x2x_(x);           %% convert real to complex x
            if nargout > 1
                [G, Gv1, Gv2, Gzr, Gzi] = obj.nodal_complex_power_balance(x_);
                Gx = [Gv1 Gv2 Gzr Gzi];
                dg = [  real(Gx);       %% P mismatch w.r.t v1, v2, zr, zi
                        imag(Gx)    ];  %% Q mismatch w.r.t v1, v2, zr, zi
            else
                G = obj.nodal_complex_power_balance(x_);
            end
            g = [ real(G);              %% active power (P) mismatch
                  imag(G) ];            %% reactive power (Q) mismatch
        end

        function d2G = opf_current_balance_hess(obj, x, lam)
            x_ = obj.x2x_(x);           %% convert real to complex x
            nlam = length(lam) / 2;
            lamIr = lam(1:nlam);
            lamIi = lam((1:nlam)+nlam);

            d2Gr = obj.nodal_complex_current_balance_hess(x_, lamIr);
            d2Gi = obj.nodal_complex_current_balance_hess(x_, lamIi);

            d2G = real(d2Gr) + imag(d2Gi);
        end

        function d2G = opf_power_balance_hess(obj, x, lam)
            x_ = obj.x2x_(x);           %% convert real to complex x
            nlam = length(lam) / 2;
            lamP = lam(1:nlam);
            lamQ = lam((1:nlam)+nlam);

            d2Gr = obj.nodal_complex_power_balance_hess(x_, lamP);
            d2Gi = obj.nodal_complex_power_balance_hess(x_, lamQ);

            d2G = real(d2Gr) + imag(d2Gi);
        end

        function add_opf_legacy_user_constraints(obj, om, mpc, mpopt)
            %% call parent
            add_opf_legacy_user_constraints@mp_aggregate(obj, om, mpc, mpopt);

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

            %% user-defined nonlinear equalities
            if nnle
                for k = 1:length(mpc.user_constraints.nle)
                    nlc = mpc.user_constraints.nle{k};
                    fcn  = eval(['@(x)' nlc{3} '(x, nlc{6}{:})']);
                    hess = eval(['@(x, lam)' nlc{4} '(x, lam, nlc{6}{:})']);
                    om.add_nln_constraint(nlc{1:2}, 1, fcn, hess, nlc{5});
                end
            end

            %% user-defined nonlinear inequalities
            if nnli
                for k = 1:length(mpc.user_constraints.nli)
                    nlc = mpc.user_constraints.nli{k};
                    fcn  = eval(['@(x)' nlc{3} '(x, nlc{6}{:})'])
                    hess = eval(['@(x, lam)' nlc{4} '(x, lam, nlc{6}{:})'])
                    om.add_nln_constraint(nlc{1:2}, 0, fcn, hess, nlc{5});
                end
            end
        end
    end     %% methods
end         %% classdef
