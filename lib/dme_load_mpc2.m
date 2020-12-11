classdef dme_load_mpc2 < dme_load & dm_format_mpc2
%DME_LOAD_MPC2  MATPOWER data model load table for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function nr = count(obj, dm)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD] = idx_bus;

            %% get bus indices
            tab = obj.get_table(dm);
            obj.bus = find(tab(:, PD) | tab(:, QD));

            %% number of loads
            nr = length(obj.bus);
            obj.nr = nr;
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elm_by_name('bus').status;  %% bus status

            %% update status of gens at isolated/offline buses
            obj.status = obj.status & bs(obj.bus);

            %% call parent to fill in on/off
            update_status@dme_load(obj, dm);
        end

        function obj = build_params(obj, dm)
            %% define named indices into data matrices
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD] = idx_bus;
            baseMVA = dm.mpc.baseMVA;

            tab = obj.get_table(dm);
            obj.Pd = tab(obj.bus(obj.on), PD) / baseMVA;
            obj.Qd = tab(obj.bus(obj.on), QD) / baseMVA;
        end

        function [s, Sd, Y] = sys_wide_zip_loads(obj, dm)
            mpc = dm.mpc;
            if isfield(mpc, 'sys_wide_zip_loads')
                pw = mpc.sys_wide_zip_loads.pw;
                qw = mpc.sys_wide_zip_loads.qw;
                if any(size(pw) ~= [1 3])
                    error('dme_load_mpc2/sys_wide_zip_loads: ''exp.sys_wide_zip_loads.pw'' must be a 1 x 3 vector');
                end
                if abs(sum(pw) - 1) > eps
                    error('dme_load_mpc2/sys_wide_zip_loads: elements of ''exp.sys_wide_zip_loads.pw'' must sum to 1');
                end
                if isempty(qw)
                    qw = pw;
                else
                    if any(size(qw) ~= [1 3])
                        error('dme_load_mpc2/sys_wide_zip_loads: ''exp.sys_wide_zip_loads.qw'' must be a 1 x 3 vector');
                    end
                    if abs(sum(qw) - 1) > eps
                        error('dme_load_mpc2/sys_wide_zip_loads: elements of ''exp.sys_wide_zip_loads.qw'' must sum to 1');
                    end
                end

                Pd = obj.Pd;
                Qd = obj.Qd;
                nd = length(Pd);

                s  = pw(1) * Pd + 1j * qw(1) * Qd;
                Sd = pw(2) * Pd + 1j * qw(2) * Qd;
                Y  = pw(3) * Pd - 1j * qw(3) * Qd;
%                i = pw(2) * Pd - 1j * qw(2) * Qd;   %% power is function of complex voltage, not voltage magnitude (as expected)
                Y = sparse(1:nd, 1:nd, Y, nd, nd);
            else
                s = []; Sd = []; Y = [];
            end
        end
    end     %% methods
end         %% classdef
