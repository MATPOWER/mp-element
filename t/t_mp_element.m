function obj = t_mp_element(quiet, out_ac)
%T_MP_ELEMENT  Tests for MP_ELEMENT.

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

t_begin(201, quiet);

define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

casefile = 't_case9_opfv2';

mpopt = mpoption('out.all', 0, 'verbose', 0);

%% DC model
t = 'dc_aggregate() : ';
dc = dc_aggregate();
t_ok(strcmp(dc.name, 'aggregate'), [t 'name']);
t_ok(strcmp(class(dc), 'dc_aggregate'), [t 'class']);
t_ok(strcmp(dc.find_model_class(), 'dc_model'), [t 'model class']);
t_ok(strcmp(dc.model_name, 'DC model'), [t 'model name']);
t_ok(strcmp(dc.model_tag, 'dc'), [t 'model tag']);
p = dc.model_params();
t_is(length(p), 3, 12, [t '# of model params']);
t_ok(strcmp(p{1}, 'B'), [t 'B']);
t_ok(strcmp(p{2}, 'K'), [t 'K']);
t_ok(strcmp(p{3}, 'p'), [t 'p']);
t_is(dc.nk, 1, 12, [t 'nk']);
t_is(dc.np, 0, 12, [t 'np']);
t_is(dc.nz, 0, 12, [t 'nz']);
t_is(dc.nv, 0, 12, [t 'nv']);
t_is(length(fieldnames(dc.set_types)), 4, 12, [t '# of set types']);
t_ok(strcmp(dc.set_types.node, 'NODES'), [t 'set_types.node']);
t_ok(strcmp(dc.set_types.state, 'STATES'), [t 'set_types.state']);
t_ok(strcmp(dc.set_types.v, 'VOLTAGE VARS (v)'), [t 'set_types.v']);
t_ok(strcmp(dc.set_types.z, 'NON-VOLTAGE VARS (z)'), [t 'set_types.z']);
t_is(length(dc.mpe_list), 0, 12, [t '# of element types']);

t = 'dc.create_model(mpc) : ';
mpc = ext2int(rundcpf(loadcase(casefile), mpopt));
t_ok(mpc.success, [t 'solved power flow']);
dc.create_model(mpc);
t_is(dc.nk, 1, 12, [t 'nk']);
t_is(dc.np, 30, 12, [t 'np']);
t_is(dc.nz, 3, 12, [t 'nz']);
t_is(dc.nv, 9, 12, [t 'nv']);
t_is(length(dc.mpe_list), 4, 12, [t '# of element types']);

mpe = dc.mpe_list;
t = 'dc_bus : '; k = 1;
t_ok(strcmp(mpe{k}.name, 'bus'), [t 'name']);
t_ok(strcmp(class(mpe{k}), 'dc_bus'), [t 'class']);
t_ok(strcmp(mpe{k}.find_model_class(), 'dc_model'), [t 'model class']);
t_ok(strcmp(mpe{k}.model_name, 'DC model'), [t 'model name']);
t_ok(strcmp(mpe{k}.model_tag, 'dc'), [t 'model tag']);
t_is(mpe{k}.nk, 9, 12, [t 'nk']);
t_is(mpe{k}.np, 0, 12, [t 'np']);
t_is(mpe{k}.nz, 0, 12, [t 'nz']);
t_ok(isempty(mpe{k}.B), [t 'B']);
t_ok(isempty(mpe{k}.K), [t 'K']);
t_ok(isempty(mpe{k}.p), [t 'p']);
% mpe{k}

t = 'dc_gen : '; k = 2;
t_ok(strcmp(mpe{k}.name, 'gen'), [t 'name']);
t_ok(strcmp(class(mpe{k}), 'dc_gen'), [t 'class']);
t_ok(strcmp(mpe{k}.find_model_class(), 'dc_model'), [t 'model class']);
t_ok(strcmp(mpe{k}.model_name, 'DC model'), [t 'model name']);
t_ok(strcmp(mpe{k}.model_tag, 'dc'), [t 'model tag']);
t_is(mpe{k}.nk, 3, 12, [t 'nk']);
t_is(mpe{k}.np, 1, 12, [t 'np']);
t_is(mpe{k}.nz, 1, 12, [t 'nz']);
t_ok(isempty(mpe{k}.B), [t 'B']);
t_is(mpe{k}.K, -speye(3), 12, [t 'K']);
t_ok(isempty(mpe{k}.p), [t 'p']);
% mpe{k}

t = 'dc_load : '; k = 3;
t_ok(strcmp(mpe{k}.name, 'load'), [t 'name']);
t_ok(strcmp(class(mpe{k}), 'dc_load'), [t 'class']);
t_ok(strcmp(mpe{k}.find_model_class(), 'dc_model'), [t 'model class']);
t_ok(strcmp(mpe{k}.model_name, 'DC model'), [t 'model name']);
t_ok(strcmp(mpe{k}.model_tag, 'dc'), [t 'model tag']);
t_is(mpe{k}.nk, 9, 12, [t 'nk']);
t_is(mpe{k}.np, 1, 12, [t 'np']);
t_is(mpe{k}.nz, 0, 12, [t 'nz']);
t_ok(isempty(mpe{k}.B), [t 'B']);
t_ok(isempty(mpe{k}.K), [t 'K']);
t_is(mpe{k}.p, [0 0 0 0 0.9 0 1 0 1.25]', 12, [t 'p']);
% mpe{k}

t = 'dc_branch : '; k = 4;
t_ok(strcmp(mpe{k}.name, 'branch'), [t 'name']);
t_ok(strcmp(class(mpe{k}), 'dc_branch'), [t 'class']);
t_ok(strcmp(mpe{k}.find_model_class(), 'dc_model'), [t 'model class']);
t_ok(strcmp(mpe{k}.model_name, 'DC model'), [t 'model name']);
t_ok(strcmp(mpe{k}.model_tag, 'dc'), [t 'model tag']);
t_is(mpe{k}.nk, 9, 12, [t 'nk']);
t_is(mpe{k}.np, 2, 12, [t 'np']);
t_is(mpe{k}.nz, 0, 12, [t 'nz']);
nl = size(mpc.branch, 1);
stat = mpc.branch(:, BR_STATUS);    %% ones at in-service branches
tap = ones(nl, 1);                  %% default tap ratio = 1
i = find(mpc.branch(:, TAP));       %% indices of non-zero tap ratios
tap(i) = mpc.branch(i, TAP);        %% assign non-zero tap ratios
b = stat ./ mpc.branch(:, BR_X);    %% series susceptance
b = b ./ tap;
Pfinj = b .* (-mpc.branch(:, SHIFT) * pi/180);
Bdc = sparse( ...
    [1:nl 1:nl nl+1:2*nl nl+1:2*nl]', ...
    [1:nl nl+1:2*nl 1:nl nl+1:2*nl]', ...
    [b; -b; -b; b], ...
    2*nl, 2*nl );
pdc = [Pfinj; -Pfinj];
t_is(mpe{k}.B, Bdc, 12, [t 'B']);
t_ok(isempty(mpe{k}.K), [t 'K']);
t_is(mpe{k}.p, pdc, 12, [t 'p']);
% mpe{k}

t = 'dc.port_inj_power(x)';
v = mpc.bus(:, VA) * pi/180;
z = mpc.gen(:, PG) / mpc.baseMVA;
x = [v;z];
P  = dc.port_inj_power(x, 1);
eP = [-0.67 -0.85 -1.63 0 0 0 0 0.9 0 1.0 0 1.25 0.67 0.289673913 -0.610326086 0.85 0.239673913 -0.760326086 -1.63 0.869673913 -0.380326086 -0.67 -0.289673913 0.610326086 -0.85 -0.239673913 0.760326086 1.63 -0.869673913 0.380326086]';
t_is(P, eP, 8, t);
C = dc.getC();
D = dc.getD();
t_is(C * P, 0, 12, [t ' : C * P == 0']);

t = 'dc.port_inj_power(A''*x, 0)';
nv = dc.nv;
nz = dc.nz;
np = dc.np;
A = [   C sparse(nv, nz);
        sparse(nz, np) D    ];
P1 = dc.port_inj_power(A'*x, 0);
t_is(P1, eP, 8, t);

t = 'dc.port_inj_power(x, 1, [3;1])';
P2  = dc.port_inj_power(x, 1, [3;1]);
t_is(P2, eP([3;1]), 8, t);

t = 'dc.mpe_by_name(''gen'') : ';
gen = dc.mpe_by_name('gen');
t_ok(strcmp(gen.name, 'gen'), [t 'name']);
t_ok(strcmp(class(gen), 'dc_gen'), [t 'class']);

t = 'gen.port_inj_power(x, 1)';
Pg = gen.port_inj_power(x, 1);
ePg = -[0.67; 0.85; 1.63];
t_is(Pg, ePg, 12, t);
 
t = 'gen.port_inj_power(x, 1, [3;1])';
Pg31 = gen.port_inj_power(x, 1, [3;1]);
t_is(Pg31, ePg([3;1]), 12, t);

t = 'gen.port_inj_power(x, 1, 2)';
Pg2 = gen.port_inj_power(x, 1, 2);
t_is(Pg2, ePg(2), 12, t);

%% AC model
t = 'acsp_aggregate() : ';
ac = acsp_aggregate();
t_ok(strcmp(ac.name, 'aggregate'), [t 'name']);
t_ok(strcmp(class(ac), 'acsp_aggregate'), [t 'class']);
t_ok(strcmp(ac.find_model_class(), 'acsp_model'), [t 'model class']);
t_ok(strcmp(ac.model_name, 'AC-power-polar model'), [t 'model name']);
t_ok(strcmp(ac.model_tag, 'acsp'), [t 'model tag']);
p = ac.model_params();
t_is(length(p), 6, 12, [t '# of model params']);
t_ok(strcmp(p{1}, 'Y'), [t 'Y']);
t_ok(strcmp(p{2}, 'L'), [t 'L']);
t_ok(strcmp(p{3}, 'M'), [t 'M']);
t_ok(strcmp(p{4}, 'N'), [t 'N']);
t_ok(strcmp(p{5}, 'i'), [t 'i']);
t_ok(strcmp(p{6}, 's'), [t 's']);
t_is(ac.nk, 1, 12, [t 'nk']);
t_is(ac.np, 0, 12, [t 'np']);
t_is(ac.nz, 0, 12, [t 'nz']);
t_is(ac.nv, 0, 12, [t 'nv']);
t_is(length(fieldnames(ac.set_types)), 6, 12, [t '# of set types']);
t_ok(strcmp(ac.set_types.node, 'NODES'), [t 'set_types.node']);
t_ok(strcmp(ac.set_types.state, 'STATES'), [t 'set_types.state']);
t_ok(strcmp(ac.set_types.va, 'VOLTAGE ANG VARS (va)'), [t 'set_types.va']);
t_ok(strcmp(ac.set_types.vm, 'VOLTAGE MAG VARS (vm)'), [t 'set_types.vm']);
t_ok(strcmp(ac.set_types.zr, 'NON-VOLTAGE VARS REAL (zr)'), [t 'set_types.zr']);
t_ok(strcmp(ac.set_types.zi, 'NON-VOLTAGE VARS IMAG (zi)'), [t 'set_types.zi']);
t_is(length(ac.mpe_list), 0, 12, [t '# of element types']);

t = 'ac.create_model(mpc) : ';
mpc = ext2int(runpf(loadcase(casefile), mpopt));
t_ok(mpc.success, [t 'solved power flow']);
ac.create_model(mpc);
t_is(ac.nk, 1, 12, [t 'nk']);
t_is(ac.np, 30, 12, [t 'np']);
t_is(ac.nz, 3, 12, [t 'nz']);
t_is(ac.nv, 18, 12, [t 'nv']);
t_is(length(ac.mpe_list), 4, 12, [t '# of element types']);

mpe = ac.mpe_list;
t = 'acsp_bus : '; k = 1;
t_ok(strcmp(mpe{k}.name, 'bus'), [t 'name']);
t_ok(strcmp(class(mpe{k}), 'acsp_bus'), [t 'class']);
t_ok(strcmp(mpe{k}.find_model_class(), 'acsp_model'), [t 'model class']);
t_ok(strcmp(mpe{k}.model_name, 'AC-power-polar model'), [t 'model name']);
t_ok(strcmp(mpe{k}.model_tag, 'acsp'), [t 'model tag']);
t_is(mpe{k}.nk, 9, 12, [t 'nk']);
t_is(mpe{k}.np, 0, 12, [t 'np']);
t_is(mpe{k}.nz, 0, 12, [t 'nz']);
t_ok(isempty(mpe{k}.Y), [t 'Y']);
t_ok(isempty(mpe{k}.L), [t 'L']);
t_ok(isempty(mpe{k}.M), [t 'M']);
t_ok(isempty(mpe{k}.N), [t 'N']);
t_ok(isempty(mpe{k}.i), [t 'i']);
t_ok(isempty(mpe{k}.s), [t 's']);
% mpe{k}

t = 'ac_gen : '; k = 2;
t_ok(strcmp(mpe{k}.name, 'gen'), [t 'name']);
t_ok(strcmp(class(mpe{k}), 'ac_gen'), [t 'class']);
t_ok(strcmp(mpe{k}.find_model_class(), 'acsp_model'), [t 'model class']);
t_ok(strcmp(mpe{k}.model_name, 'AC-power-polar model'), [t 'model name']);
t_ok(strcmp(mpe{k}.model_tag, 'acsp'), [t 'model tag']);
t_is(mpe{k}.nk, 3, 12, [t 'nk']);
t_is(mpe{k}.np, 1, 12, [t 'np']);
t_is(mpe{k}.nz, 1, 12, [t 'nz']);
t_ok(isempty(mpe{k}.Y), [t 'Y']);
t_ok(isempty(mpe{k}.L), [t 'L']);
t_ok(isempty(mpe{k}.M), [t 'M']);
t_is(mpe{k}.N, -speye(3), 12, [t 'N']);
t_ok(isempty(mpe{k}.i), [t 'i']);
t_ok(isempty(mpe{k}.s), [t 's']);
% mpe{k}

t = 'ac_load : '; k = 3;
t_ok(strcmp(mpe{k}.name, 'load'), [t 'name']);
t_ok(strcmp(class(mpe{k}), 'ac_load'), [t 'class']);
t_ok(strcmp(mpe{k}.find_model_class(), 'acsp_model'), [t 'model class']);
t_ok(strcmp(mpe{k}.model_name, 'AC-power-polar model'), [t 'model name']);
t_ok(strcmp(mpe{k}.model_tag, 'acsp'), [t 'model tag']);
t_is(mpe{k}.nk, 9, 12, [t 'nk']);
t_is(mpe{k}.np, 1, 12, [t 'np']);
t_is(mpe{k}.nz, 0, 12, [t 'nz']);
t_ok(isempty(mpe{k}.Y), [t 'Y']);
t_ok(isempty(mpe{k}.L), [t 'L']);
t_ok(isempty(mpe{k}.M), [t 'M']);
t_ok(isempty(mpe{k}.N), [t 'N']);
t_ok(isempty(mpe{k}.i), [t 'i']);
t_is(mpe{k}.s, [0 0 0 0 0.9 0 1 0 1.25]' + 1j*[0 0 0 0 0.3 0 0.35 0 0.5]', 12, [t 's']);
% mpe{k}

t = 'ac_branch : '; k = 4;
t_ok(strcmp(mpe{k}.name, 'branch'), [t 'name']);
t_ok(strcmp(class(mpe{k}), 'ac_branch'), [t 'class']);
t_ok(strcmp(mpe{k}.find_model_class(), 'acsp_model'), [t 'model class']);
t_ok(strcmp(mpe{k}.model_name, 'AC-power-polar model'), [t 'model name']);
t_ok(strcmp(mpe{k}.model_tag, 'acsp'), [t 'model tag']);
t_is(mpe{k}.nk, 9, 12, [t 'nk']);
t_is(mpe{k}.np, 2, 12, [t 'np']);
t_is(mpe{k}.nz, 0, 12, [t 'nz']);
nl = size(mpc.branch, 1);
stat = mpc.branch(:, BR_STATUS);    %% ones at in-service branches
tap = ones(nl, 1);                  %% default tap ratio = 1
i = find(mpc.branch(:, TAP));       %% indices of non-zero tap ratios
tap(i) = mpc.branch(i, TAP);        %% assign non-zero tap ratios
tap = tap .* exp(1j*pi/180 * mpc.branch(:, SHIFT)); %% add phase shifters
Ys = stat ./ (mpc.branch(:, BR_R) + 1j * mpc.branch(:, BR_X));  %% series admittance
Bc = stat .* mpc.branch(:, BR_B);   %% line charging susceptance
Ytt = Ys + 1j*Bc/2;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;
eY = sparse( ...
    [1:nl 1:nl nl+1:2*nl nl+1:2*nl]', ...
    [1:nl nl+1:2*nl 1:nl nl+1:2*nl]', ...
    [Yff; Yft; Ytf; Ytt], 2*nl, 2*nl );
t_is(mpe{k}.Y, eY, 12, [t 'Y']);
t_ok(isempty(mpe{k}.L), [t 'L']);
t_ok(isempty(mpe{k}.M), [t 'M']);
t_ok(isempty(mpe{k}.N), [t 'N']);
t_ok(isempty(mpe{k}.i), [t 'i']);
t_ok(isempty(mpe{k}.s), [t 's']);
% mpe{k}

t = 'S = ac.port_inj_power(x)';
v = mpc.bus(:, VM) .* exp(1j * mpc.bus(:, VA) * pi/180);
z = (mpc.gen(:, PG) + 1j * mpc.gen(:, QG)) / mpc.baseMVA;
x = [v;z];
S  = ac.port_inj_power(x);
eS = [-0.7195470 -0.85 -1.63 0 0 0 0 0.9 0 1.0 0 1.25 0.7195470 0.3072828 -0.5944531 0.85 0.2410613 -0.7598935 -1.63 0.8650443 -0.4096011 -0.7195470 -0.3055468 0.6089386 -0.85 -0.2401064 0.7649556 1.63 -0.8403988 0.4122642]' ...
  + 1j * [-0.2406895 0.0364902 -0.1446011 0 0 0 0 0.3 0 0.35 0 0.5 0.2406895 -0.0058585 -0.1631204 -0.0364902 0.0453679 -0.1059923 0.0227618 -0.0253242 -0.3571801 -0.2075304 -0.1368795 -0.1242746 0.0789067 -0.2440076 0.0025623 0.1446011 -0.1428198 0.2133889]';
t_is(S, eS, 6, t);
C = ac.getC();
D = ac.getD();
t_is(C * S, 0, 10, [t ' : C * S == 0']);

t = 'S = ac.port_inj_power(A''*x, 0)';
nv = ac.nv;
nz = ac.nz;
np = ac.np;
A = [   C sparse(nv/2, nz);
        sparse(nz, np) D    ];
S1 = ac.port_inj_power(A'*x, 0);
t_is(S1, eS, 6, t);

t = 'ac.port_inj_power(x, 1, [3;1])';
S2  = ac.port_inj_power(x, 1, [3;1]);
t_is(S2, eS([3;1]), 6, t);

t = '[S, Sva, Svm, Szr, Szi] = ac.port_inj_power(x, 1) : ';
[S, Sva, Svm, Szr, Szi] = ac.port_inj_power(x, 1);
t_is(S, eS, 6, [t 'S']);
t_is(size(Sva), [30 9], 12, [t 'size(Sva)']);
t_is(size(Svm), [30 9], 12, [t 'size(Svm)']);
t_is(size(Szr), [30 3], 12, [t 'size(Szr)']);
t_is(size(Szi), [30 3], 12, [t 'size(Szi)']);

t = '[S, Sva, Svm, Szr, Szi] = ac.port_inj_power(x, 1, [3;2;1]) : ';
[S1, Sva1, Svm1, Szr1, Szi1] = ac.port_inj_power(x, 1, [3;2;1]);
t_is(S1, eS([3;2;1]), 6, [t 'S']);
t_is(Sva1, Sva([3;2;1], :), 12, [t 'Sva']);
t_is(Svm1, Svm([3;2;1], :), 12, [t 'Svm']);
t_is(Szr1, Szr([3;2;1], :), 12, [t 'Szr']);
t_is(Szi1, Szi([3;2;1], :), 12, [t 'Szi']);

t = 'ac.mpe_by_name(''gen'') : ';
gen = ac.mpe_by_name('gen');
t_ok(strcmp(gen.name, 'gen'), [t 'name']);
t_ok(strcmp(class(gen), 'ac_gen'), [t 'class']);

t = 'gen.port_inj_power(x, 1)';
Sg = gen.port_inj_power(x, 1);
eSg = -[0.7195470; 0.85; 1.63] + 1j * [-0.2406895; 0.0364902; -0.1446011];
t_is(Sg, eSg, 6, t);
 
t = 'gen.port_inj_power(x, 1, [3;1])';
Sg31 = gen.port_inj_power(x, 1, [3;1]);
t_is(Sg31, eSg([3;1]), 6, t);

t = 'gen.port_inj_power(x, 1, 2)';
Sg2 = gen.port_inj_power(x, 1, 2);
t_is(Sg2, eSg(2), 6, t);

t = 'H = ac.port_inj_power_hess(x, lam, 1) : ';
lam = [1:np]' / 100;
H = ac.port_inj_power_hess(x, lam, 1);
t_is(size(H), [24 24], 12, t);

t = 'H = ac.port_inj_power_hess(x, lam(1:3), 1, [3;2;1]) : ';
H = ac.port_inj_power_hess(x, lam(1:3), 1, [3;2;1]);
t_is(size(H), [24 24], 12, t);

t = 'I = ac.port_inj_current(x)';
I  = ac.port_inj_current(x);
eI = [-0.7195470 -0.8440196 -1.6311322 0 0 0 0 0.8988176 0 1.0183565 0 1.2619557 0.7195470 0.3113025 -0.5961880 0.8440196 0.2416340 -0.7720819 -1.6311322 0.8647643 -0.3982054 -0.7195470 -0.3026295 0.6023856 -0.8440196 -0.2462745 0.7663679 1.6311322 -0.8637503 0.4082444]' ...
  + 1j * [0.2406895 -0.1070623 -0.1312139 0 0 0 0 -0.3714244 0 -0.3440707 0 -0.6196286 -0.2406895 -0.0071427 0.2095040 0.1070623 -0.0371169 0.0991665 -0.1312139 0.0829015 0.4043547 0.2406895 0.1619204 0.1441792 -0.1070623 0.2449042 0.0483124 0.1312139 0.2152738 -0.2335468]';
t_is(I, eI, 6, t);
t_is(C * I, 0, 10, [t ' : C * I == 0']);

t = 'I = ac.port_inj_current(A''*x, 0)';
I1 = ac.port_inj_current(A'*x, 0);
t_is(I1, eI, 6, t);

t = '[I, Iva, Ivm, Izr, Izi] = ac.port_inj_current(x, 1) : ';
[I, Iva, Ivm, Izr, Izi] = ac.port_inj_current(x, 1);
t_is(I, eI, 6, [t 'I']);
t_is(size(Iva), [30 9], 12, [t 'size(Iva)']);
t_is(size(Ivm), [30 9], 12, [t 'size(Ivm)']);
t_is(size(Izr), [30 3], 12, [t 'size(Izr)']);
t_is(size(Izi), [30 3], 12, [t 'size(Izi)']);

t = '[I, Iva, Ivm, Izr, Izi] = ac.port_inj_current(x, 1, [3;2;1]) : ';
[I1, Iva1, Ivm1, Izr1, Izi1] = ac.port_inj_current(x, 1, [3;2;1]);
t_is(I1, eI([3;2;1]), 6, [t 'I']);
t_is(Iva1, Iva([3;2;1], :), 12, [t 'Iva']);
t_is(Ivm1, Ivm([3;2;1], :), 12, [t 'Ivm']);
t_is(Izr1, Izr([3;2;1], :), 12, [t 'Izr']);
t_is(Izi1, Izi([3;2;1], :), 12, [t 'Izi']);

%% AC Newton power flow
t = '[x, success, i] = ac.solve_power_flow(mpc, mpopt) : ';
mpc = ext2int(loadcase(casefile));
mpc = rmfield(mpc, 'order');
ac = acsp_aggregate().create_model(mpc);
[x, success, i] = ac.solve_power_flow(mpc, mpopt);
ex = [0.168751367 0.083270936 -0.042003860 -0.070114489 0.033608089 0.010847998 0.066307156 -0.075920663 0.987006852 0.975472177 1.003375436 0.985644881 0.996185245 0.957621040]';
t_is(x, ex, 8, [t 'x']);
t_is(success, 1, 12, [t 'success']);
t_is(i, 4, 12, [t 'i']);

% disp(dc)
% disp(ac)

t_end;

if nargout
    if nargin < 2
        out_ac = 0;
    end
    if out_ac
        obj = ac;
    else
        obj = dc;
    end
end

return;
