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
                obj.create_model();
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
            if nargout > 1
                dm1 = obj.copy();
                dm1.mpc.bus(:, BS) = 0;         %% zero out shunts at buses
                dm2 = dm1.copy();
                dm1.mpc.branch(:, BR_B) = 0;    %% zero out line charging shunts
                dm1.mpc.branch(:, TAP) = 1;     %% cancel out taps
                if strcmp(alg, 'FDXB')          %% if XB method
                    dm1.mpc.branch(:, BR_R) = 0;%% zero out line resistance
                end
            else
                dm2 = obj.copy();
            end

            %% modify data model to form Bpp (B double prime)
            dm2.mpc.branch(:, SHIFT) = 0;   %% zero out phase shifters
            if strcmp(alg, 'FDBX')          %% if BX method
                dm2.mpc.branch(:, BR_R) = 0;%% zero out line resistance
            end

            if nargout == 1
                dm1 = dm2;
            end
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
    end     %% methods
end         %% classdef
