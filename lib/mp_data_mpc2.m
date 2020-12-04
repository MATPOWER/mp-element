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
        function obj = mp_data_mpc2(m)
            %% call parent constructor
            obj@mp_data();

            if nargin
                %% define constants
                [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;
                [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS] = idx_gen;
                [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                    TAP, SHIFT, BR_STATUS] = idx_brch;

                %% assign table names, id and status columns
                obj.tab = struct( ...
                        'name',     {'bus', 'gen', 'branch'}, ...
                        'id_col',   {BUS_I, 0, 0}, ...
                        'st_col',   {BUS_TYPE, GEN_STATUS, BR_STATUS} ...
                    );

                %% load case and create mappings
                obj.mpc = loadcase(m);
                obj.create_model();
            end
        end

        function [ids, N, M] = get_ids(obj, tab_name, id_col)
            data = obj.mpc.(tab_name);
            if id_col && isfield(obj.mpc, tab_name)
                ids = data(:, id_col);
            else
                ids = [];
            end
            N = size(data, 1);  %% number of rows (length of di2id)
            M = max(ids);       %% max index (length of id2di)
        end

        function stats = get_status(obj, tab_name, st_col)
            switch tab_name
                case 'bus'
                    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE] = idx_bus;

                    %% check that all buses have a valid BUS_TYPE
                    bt = obj.mpc.bus(:, st_col);
                    err = find(~(bt == PQ | bt == PV | bt == REF | bt == NONE));
                    if ~isempty(err)
                        error('ext2int: bus %d has an invalid BUS_TYPE', err);
                    end
                    stats = (bt ~= NONE);       %% bus status
                otherwise
                    if st_col
                        stats = obj.mpc.(tab_name)(:, st_col) > 0;
                    else
                        stats = [];
                    end
            end
        end

        function obj = update_status(obj)
            %% define constants
            [GEN_BUS] = idx_gen;
            [F_BUS, T_BUS] = idx_brch;

            bm = obj.map.bus;               %% bus mapping info
            bs = bm.status;                 %% bus status
            b2i = bm.id2di;                 %% bus num to idx mapping

            %% update status of gen at isolated/offline bus
            gs = obj.map.gen.status;        %% gen status
            gs = gs & bs(b2i(obj.mpc.gen(:, GEN_BUS)));
            obj.map.gen.status = gs;

            %% update status of branch connected to isolated/offline bus
            brs = obj.map.branch.status;    %% branch status
            brs = brs & bs(b2i(obj.mpc.branch(:, F_BUS))) & ...
                        bs(b2i(obj.mpc.branch(:, T_BUS)));
            obj.map.branch.status = brs;

            %% call parent to fill in on/off
            update_status@mp_data(obj);
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
