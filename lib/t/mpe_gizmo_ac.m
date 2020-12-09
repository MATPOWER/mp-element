classdef mpe_gizmo_ac < mpe_gizmo% & mp_model_ac

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'gizmo';
%     end
    
    methods
        function obj = add_zvars(obj, nm, dm, idx)
            mpc = dm.mpc;
            nk = obj.nk;
            switch idx{:}
                case 1
                    Zmax = ones(nk, 1);
                    Zr   = mpc.gizmo(:, 20);
                    Zi   = mpc.gizmo(:, 21);
                case 2
                    Zmax = 2 * ones(nk, 1);
                    Zr   = mpc.gizmo(:, 22);
                    Zi   = mpc.gizmo(:, 23);
            end
            vname_r = sprintf('Zr%d_gizmo', idx{:});
            vname_i = sprintf('Zi%d_gizmo', idx{:});
            nm.add_var('zr', vname_r, nk, Zr, -Zmax, Zmax);
            nm.add_var('zi', vname_i, nk, Zi, -Zmax, Zmax);
        end

        function obj = build_params(obj, nm, dm)
            build_params@mpe_gizmo(obj, nm, dm);   %% call parent
            nk = obj.nk;
            mpc = dm.mpc;

            %% collect parameters from mpc
            y1 = mpc.gizmo(:,  4) + 1j * mpc.gizmo(:,  5);
            y2 = mpc.gizmo(:,  6) + 1j * mpc.gizmo(:,  7);
            ll = mpc.gizmo(:,  8) + 1j * mpc.gizmo(:,  9);
            ii = mpc.gizmo(:, 10) + 1j * mpc.gizmo(:, 11);
            m1 = mpc.gizmo(:, 12) + 1j * mpc.gizmo(:, 13);
            m2 = mpc.gizmo(:, 14) + 1j * mpc.gizmo(:, 15);
            nn = mpc.gizmo(:, 16) + 1j * mpc.gizmo(:, 17);
            ss = mpc.gizmo(:, 18) + 1j * mpc.gizmo(:, 19);
            zz = zeros(nk, 1);

            %% construct model parameters
            j1 = (1:nk);
            j2 = nk+j1;
            j3 = nk+j2;
            obj.Y = sparse( ...
                [j1 j1 j1 j2 j2 j2 j3 j3 j3]', ...
                [j1 j2 j3 j1 j2 j3 j1 j2 j3]', ...
                [y1; zz; -y1; zz; y2; zz; -y1; zz; y1], 3*nk, 3*nk );
            obj.L = sparse( ...
                [j1 j1 j2 j2 j3 j3 ]', ...
                [j1 j2 j1 j2 j1 j2 ]', ...
                [zz; ll; zz; -ll; zz; zz], 3*nk, 2*nk );
            obj.i = [-ii; ii; zz];
            obj.M = sparse( ...
                [j1 j1 j1 j2 j2 j2 j3 j3 j3]', ...
                [j1 j2 j3 j1 j2 j3 j1 j2 j3]', ...
                [m1; -m1; zz; -m1; m1; zz; zz; zz; m2], 3*nk, 3*nk );
            obj.N = sparse( ...
                [j1 j1 j2 j2 j3 j3 ]', ...
                [j1 j2 j1 j2 j1 j2 ]', ...
                [zz; zz; nn; zz; -nn; zz], 3*nk, 2*nk );
            obj.s = [zz; -ss; ss];
        end
    end     %% methods
end         %% classdef
