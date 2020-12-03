function obj = t_mp_data(quiet, out_ac)
%T_MP_DATA  Tests for MP_DATA.

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

% define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

%% define constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

casefile = 't_case_ext';
mpc = loadcase(casefile);
nb = size(mpc.bus, 1);
id2di = zeros(2, 1);    %% initialize as col vector
id2di(mpc.bus(:, BUS_I)) = (1:nb);

tests = {
    {'filename', casefile},
    {'mpc', mpc}
};
nt = length(tests);

t_begin(32*nt, quiet);

for k = 1:nt
    t = sprintf('mp_data_mpc2(%s) : ', tests{k}{1});
    dm = mp_data_mpc2(tests{k}{2});
    t_ok(strcmp(dm.mpc.version, '2'), [t 'mpc version']);
    t_ok(isfield(dm.tab, 'name'), [t 'dm.tab.name exists']);
    t_ok(isfield(dm.tab, 'id_col'), [t 'dm.tab.id_col exists']);
    t_ok(isfield(dm.tab, 'st_col'), [t 'dm.tab.st_col exists']);
    t_is(length(dm.tab), 3, 12, [t 'length(dm.tab)']);
    t_ok(isfield(dm.map, 'bus'), [t 'dm.map.bus exists']);
    t_ok(isfield(dm.map, 'gen'), [t 'dm.map.gen exists']);
    t_ok(isfield(dm.map, 'branch'), [t 'dm.map.branch exists']);

    t_is(dm.map.bus.k, 1, 12, [t 'bus.k']);
    t_is(dm.map.bus.N, 10, 12, [t 'bus.N']);
    t_is(dm.map.bus.n, 9, 12, [t 'bus.n']);
    t_is(dm.map.bus.di2id, mpc.bus(:, BUS_I), 12, [t 'bus.di2id']);
    t_is(dm.map.bus.id2di, id2di, 12, [t 'bus.id2di']);
    t_is(dm.map.bus.status, [1;1;1;1;1;0;1;1;1;1], 12, [t 'bus.status']);
    t_is(dm.map.bus.on, [1;2;3;4;5;7;8;9;10], 12, [t 'bus.on']);
    t_is(dm.map.bus.off, 6, 12, [t 'bus.off']);

    t_is(dm.map.gen.k, 2, 12, [t 'gen.k']);
    t_is(dm.map.gen.N, 4, 12, [t 'gen.N']);
    t_is(dm.map.gen.n, 3, 12, [t 'gen.n']);
    t_ok(isempty(dm.map.gen.di2id), [t 'gen.di2id']);
    t_ok(isempty(dm.map.gen.id2di), [t 'gen.id2di']);
    t_is(dm.map.gen.status, [1;1;0;1], 12, [t 'gen.status']);
    t_is(dm.map.gen.on, [1;2;4], 12, [t 'gen.on']);
    t_is(dm.map.gen.off, 3, 12, [t 'gen.off']);

    t_is(dm.map.branch.k, 3, 12, [t 'branch.k']);
    t_is(dm.map.branch.N, 10, 12, [t 'branch.N']);
    t_is(dm.map.branch.n, 9, 12, [t 'branch.n']);
    t_ok(isempty(dm.map.branch.di2id), [t 'branch.di2id']);
    t_ok(isempty(dm.map.branch.id2di), [t 'branch.id2di']);
    t_is(dm.map.branch.status, [1;1;1;1;1;1;0;1;1;1], 12, [t 'branch.status']);
    t_is(dm.map.branch.on, [1;2;3;4;5;6;8;9;10], 12, [t 'branch.on']);
    t_is(dm.map.branch.off, 7, 12, [t 'branch.off']);
end

t_end;
