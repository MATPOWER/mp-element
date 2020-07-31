classdef ac_load < mp_load% & ac_model

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'load';
%     end

    methods
        function k = load_bus(obj, mpc)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            k = find(mpc.bus(:, PD) | mpc.bus(:, QD));
        end

        function obj = build_params(obj, nm, mpc)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
                VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

            build_params@mp_load(obj, nm, mpc);    %% call parent

            obj.s = (mpc.bus(obj.busidx, PD) + ...
                1j * mpc.bus(obj.busidx, QD)) / mpc.baseMVA;    %% vector of complex power demand

            %%-----  HACK ALERT  -----
            %% This is a hack to deal with experimental
            %% mpopt.exp.sys_wide_zip_loads.pw/qw.
            if isfield(mpc, 'sys_wide_zip_loads')
                pw = mpc.sys_wide_zip_loads.pw;
                qw = mpc.sys_wide_zip_loads.qw;
                if any(size(pw) ~= [1 3])
                    error('''exp.sys_wide_zip_loads.pw'' must be a 1 x 3 vector');
                end
                if abs(sum(pw) - 1) > eps
                    error('elements of ''exp.sys_wide_zip_loads.pw'' must sum to 1');
                end
                if isempty(qw)
                    qw = pw;
                else
                    if any(size(qw) ~= [1 3])
                        error('''exp.sys_wide_zip_loads.qw'' must be a 1 x 3 vector');
                    end
                    if abs(sum(qw) - 1) > eps
                        error('elements of ''exp.sys_wide_zip_loads.qw'' must sum to 1');
                    end
                end

                Pd = mpc.bus(obj.busidx, PD) / mpc.baseMVA;
                Qd = mpc.bus(obj.busidx, QD) / mpc.baseMVA;
                nd = length(Pd);

                obj.s = pw(1) * Pd + 1j * qw(1) * Qd;
                Sd    = pw(2) * Pd + 1j * qw(2) * Qd;
                Y     = pw(3) * Pd - 1j * qw(3) * Qd;
%                obj.i = pw(2) * Pd - 1j * qw(2) * Qd;   %% power is function of complex voltage, not voltage magnitude (as expected)
                obj.Y = sparse(1:nd, 1:nd, Y, nd, nd);
                obj.snln = @(x_, sysx, idx)port_inj_power_nln(obj, Sd, x_, sysx, idx);
            end
            %%-----  end of HACK  -----

            function [S, Sv1, Sv2, Szr, Szi] = port_inj_power_nln(obj, Sd, x_, sysx, idx)
                if nargin < 5
                    idx = [];
                    if nargin < 4
                        sysx = 1;
                    end
                end

                [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);
                if isempty(idx)
                    Sdi = Sd;
                else
                    Sdi = Sd(idx);
                end
                S = abs(vi_) .* Sdi;

                if nargout > 1
                    nv = length(v_);
                    nz = length(z_);
                    ni = length(S);
                    dSdi = spdiags(Sdi, 0, ni, ni);
                    if isempty(idx)
                        idx = (1:ni);
                    end
                    Sv1 = sparse(ni, nv);
                    Sv2 = sparse(1:ni, idx, Sdi, ni, nv);
                    if nargout > 3
                        Szr = sparse(ni, nz);
                        Szi = Szr;
                    end
                    if sysx
                        Ct = obj.C';
                        Sv1 = Sv1 * Ct;
                        Sv2 = Sv2 * Ct;
                        if nargout > 3  %% Szr, Szi are empty, but num of rows is needed
                            Dt = obj.D';
                            Szr = Szr * Dt;
                            Szi = Szi * Dt;
                        end
                    end
                end
            end
        end
    end     %% methods
end         %% classdef
