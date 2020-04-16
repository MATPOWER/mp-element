classdef dc_aggregate < mp_aggregate & dc_model

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        va = [];
        z  = [];
    end
    
    methods
        %% constructor
        function obj = dc_aggregate()
            obj@mp_aggregate();
            obj.element_classes = { @dc_bus, @dc_gen, @dc_load, @dc_branch, @dc_shunt };
            if isempty(obj.node)    %% skip if constructed from existing object
                obj.init_set_types();   %% should be called in mp_idx_manager
                                        %% constructor, if not for:
            end                         %% https://savannah.gnu.org/bugs/?52614
        end

        function obj = def_set_types(obj)
            def_set_types@mp_aggregate(obj);        %% call parent first
            obj.set_types.va = 'VOLTAGE VARS (va)';
            obj.set_types.z  = 'NON-VOLTAGE VARS (z)';
        end

        function obj = build_params(obj, asm, mpc)
            %% call parent to build individual element parameters
            build_params@mp_aggregate(obj, asm, mpc);

            %% aggregate parameters from individual elements
            obj.B = obj.stack_matrix_params('B', 1);
            obj.K = obj.stack_matrix_params('K', 0);
            obj.p = obj.stack_vector_params('p');
        end

        function add_opf_node_balance_constraints(obj, om)
            [B, K, p] = obj.get_params();

            %% power balance constraints
            C = obj.C;
            Amis = C * [B*C' K*obj.D'];
            bmis = -C * p;
            om.add_lin_constraint('Pmis', Amis, bmis, bmis, ...
                                {obj.va.order(:).name obj.z.order(:).name});

            %% user data
            branch_mpe = obj.mpe_by_name('branch');
            [Bbr, pbr] = branch_mpe.get_params(1:branch_mpe.nk, {'B', 'p'});
            om.userdata.Bf = Bbr * branch_mpe.C';
            om.userdata.Pfinj = pbr;
        end

        function names = opf_legacy_user_var_names(obj)
            names = {'Va', 'Pg'};
        end

        function [va, success, i, data] = solve_power_flow(obj, mpc, mpopt)
            %% MATPOWER options
            if nargin < 3
                mpopt = mpoption;
            end

            if mpopt.verbose, fprintf('-----  solve_power_flow()  -----\n'); end

            %% get bus index lists of each type of bus
            [ref, pv, pq] = bustypes(mpc.bus, mpc.gen);
            npv = length(pv);
            npq = length(pq);

            %% create x0 for Newton power flow
            va = obj.params_var('va');
            z = obj.params_var('z');

            %% constant
            va_threshold = 1e5;     %% arbitrary threshold on |va| for declaring failure
            i = 1;                  %% not iterative

            %% set up to trap non-singular matrix warnings
            [lastmsg, lastid] = lastwarn;
            lastwarn('');

            %% get parameters
            [B, K, p] = obj.get_params();
            BB = obj.C * B * obj.C';
            pbus = -(obj.C * K * obj.D' * z + obj.C * p);

            %% update angles for non-reference buses
            va([pv; pq]) = BB([pv; pq], [pv; pq]) \ ...
                            (pbus([pv; pq]) - BB([pv; pq], ref) * va(ref));

            [msg, id] = lastwarn;
            %% Octave is not consistent in assigning proper warning id, so we'll just
            %% check for presence of *any* warning
            success = 1;    %% successful by default
            if ~isempty(msg) || max(abs(va)) > va_threshold
                success = 0;
            end

            %% restore warning state
            lastwarn(lastmsg, lastid);

            branch_mpe = obj.mpe_by_name('branch');
            [Bf, pf] = branch_mpe.get_params(1:branch_mpe.nk, {'B', 'p'});
            data = struct( ...
                    'B', BB, ...
                    'Bf', Bf * branch_mpe.C', ...
                    'Pbus', pbus, ...
                    'Pfinj', pf  ...
                );
        end
    end     %% methods
end         %% classdef
