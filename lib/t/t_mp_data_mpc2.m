function obj = t_mp_data_mpc2(quiet)
%T_MP_DATA_MPC2  Tests for MP_DATA_MPC2.

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

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

casefile = 't_case_ext';
mpc = loadcase(casefile);
dm0 = mp_data_mpc2().build(mpc);
nb = size(mpc.bus, 1);
ID2i = zeros(2, 1);     %% initialize as col vector
ID2i(mpc.bus(:, BUS_I)) = (1:nb);

tests = {
    {'filename', casefile},
    {'mpc', mpc},
    {'dm', dm0}
};
nt = length(tests);

t_begin(66*nt, quiet);

for k = 1:nt
    if isa(tests{k}{2}, 'mp_data')
        t = sprintf('%s.copy() : ', tests{k}{1});
        dm = tests{k}{2}.copy();
        tests{k}{2}.mpc = [];
        tests{k}{2}.elm_list{1}.ID = [];
    else
        t = sprintf('mp_data_mpc2().build(%s) : ', tests{k}{1});
        dm = mp_data_mpc2().build(tests{k}{2});
    end
    t_ok(strcmp(dm.mpc.version, '2'), [t 'mpc version']);
    t_ok(iscell(dm.element_classes), [t 'iscell(dm.element_classes)']);
    t_is(length(dm.element_classes), 5, 12, [t 'length(dm.element_classes)']);
    t_ok(iscell(dm.elm_list), [t 'iscell(dm.elm_list)']);
    t_is(length(dm.elm_list), 5, 12, [t 'length(dm.elm_list)']);
    t_ok(isstruct(dm.elm_map), [t 'isstruct(dm.elm_map)']);
    t_ok(isfield(dm.elm_map, 'bus'), [t 'dm.elm_map.bus exists']);
    t_ok(isfield(dm.elm_map, 'gen'), [t 'dm.elm_map.gen exists']);
    t_ok(isfield(dm.elm_map, 'load'), [t 'dm.elm_map.load exists']);
    t_ok(isfield(dm.elm_map, 'branch'), [t 'dm.elm_map.branch exists']);
    t_ok(isfield(dm.elm_map, 'shunt'), [t 'dm.elm_map.shunt exists']);
    t_is(dm.elm_map.bus, 1, 12, [t 'dm.elm_map.bus']);
    t_is(dm.elm_map.gen, 2, 12, [t 'dm.elm_map.gen']);
    t_is(dm.elm_map.load, 3, 12, [t 'dm.elm_map.load']);
    t_is(dm.elm_map.branch, 4, 12, [t 'dm.elm_map.branch']);
    t_is(dm.elm_map.shunt, 5, 12, [t 'dm.elm_map.shunt']);

    bus = dm.elm_by_name('bus');
    t_ok(isa(bus, 'dme_bus_mpc2'), [t 'bus class']);
    t_ok(isa(bus, 'dm_element'), [t 'bus isa dm_element']);
    t_ok(strcmp(bus.name, 'bus'), [t 'bus.name']);
    t_is(bus.nr, 10, 12, [t 'bus.nr']);
    t_is(bus.n, 9, 12, [t 'bus.n']);
    t_is(bus.ID, mpc.bus(:, BUS_I), 12, [t 'bus.ID']);
    t_is(bus.ID2i, ID2i, 12, [t 'bus.ID2i']);
    t_is(bus.status, [1;1;1;1;1;0;1;1;1;1], 12, [t 'bus.status']);
    t_is(bus.on, [1;2;3;4;5;7;8;9;10], 12, [t 'bus.on']);
    t_is(bus.off, 6, 12, [t 'bus.off']);

    gen = dm.elm_by_name('gen');
    t_ok(isa(gen, 'dme_gen_mpc2'), [t 'gen class']);
    t_ok(isa(gen, 'dm_element'), [t 'gen isa dm_element']);
    t_ok(strcmp(gen.name, 'gen'), [t 'gen.name']);
    t_is(gen.nr, 4, 12, [t 'gen.nr']);
    t_is(gen.n, 3, 12, [t 'gen.n']);
    t_ok(isempty(gen.ID), [t 'gen.ID']);
    t_ok(isempty(gen.ID2i), [t 'gen.ID2i']);
    t_is(gen.status, [1;1;0;1], 12, [t 'gen.status']);
    t_is(gen.on, [1;2;4], 12, [t 'gen.on']);
    t_is(gen.off, 3, 12, [t 'gen.off']);

    ld = dm.elm_by_name('load');
    t_ok(isa(ld, 'dme_load_mpc2'), [t 'load class']);
    t_ok(isa(ld, 'dm_element'), [t 'ld isa dm_element']);
    t_ok(strcmp(ld.name, 'load'), [t 'ld.name']);
    t_is(ld.nr, 3, 12, [t 'ld.nr']);
    t_is(ld.n, 3, 12, [t 'ld.n']);
    t_ok(isempty(ld.ID), [t 'ld.ID']);
    t_ok(isempty(ld.ID2i), [t 'ld.ID2i']);
    t_is(ld.status, [1;1;1], 12, [t 'ld.status']);
    t_is(ld.on, [1;2;3], 12, [t 'ld.on']);
    t_ok(isempty(ld.off), [t 'ld.off']);

    branch = dm.elm_by_name('branch');
    t_ok(isa(branch, 'dme_branch_mpc2'), [t 'branch class']);
    t_ok(isa(branch, 'dm_element'), [t 'branch isa dm_element']);
    t_ok(strcmp(branch.name, 'branch'), [t 'branch.name']);
    t_is(branch.nr, 10, 12, [t 'branch.nr']);
    t_is(branch.n, 9, 12, [t 'branch.n']);
    t_ok(isempty(branch.ID), [t 'branch.ID']);
    t_ok(isempty(branch.ID2i), [t 'branch.ID2i']);
    t_is(branch.status, [1;1;1;1;1;1;0;1;1;1], 12, [t 'branch.status']);
    t_is(branch.on, [1;2;3;4;5;6;8;9;10], 12, [t 'branch.on']);
    t_is(branch.off, 7, 12, [t 'branch.off']);

    shunt = dm.elm_by_name('shunt');
    t_ok(isa(shunt, 'dme_shunt_mpc2'), [t 'shunt class']);
    t_ok(isa(shunt, 'dm_element'), [t 'shunt isa dm_element']);
    t_ok(strcmp(shunt.name, 'shunt'), [t 'shunt.name']);
    t_is(shunt.nr, 2, 12, [t 'shunt.nr']);
    t_is(shunt.n, 2, 12, [t 'shunt.n']);
    t_ok(isempty(shunt.ID), [t 'shunt.ID']);
    t_ok(isempty(shunt.ID2i), [t 'shunt.ID2i']);
    t_is(shunt.status, [1;1], 12, [t 'shunt.status']);
    t_is(shunt.on, [1;2], 12, [t 'shunt.on']);
    t_ok(isempty(shunt.off), [t 'shunt.off']);
end

t_end;
