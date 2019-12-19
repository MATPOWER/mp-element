function obj = t_acp_port_inj_power(quiet)
%T_ACP_PORT_INJ_POWER  Tests of port_inj_power() derivatives using MP_ELEMENT.

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

t_begin(87, quiet);

define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

casefile = 't_case9_opfv2';

mpopt = mpoption('out.all', 0, 'verbose', 0);

%% create aggregate system object
mpc = ext2int(loadcase(casefile));
mpc = rmfield(mpc, 'order');
ac = acsp_aggregate().create_model(mpc);
C = ac.getC();
D = ac.getD();
np = ac.np;
nv = ac.nv/2;
nz = ac.nz;
A = [   C sparse(nv, nz);
        sparse(nz, np) D    ];
A2 = [A sparse(nv+nz, np+nz); sparse(nv+nz, np+nz) A];

%% other parameters
dx = 1e-8;
idx = randperm(np, fix(0.67*np))';
lam = (1.5*rand(np, 1) + 0.5); k = randperm(np, fix(np/2)); lam(k) = -lam(k);
e0 = zeros(np, 1);
e1 = ones(np, 1);

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(mpc.bus, mpc.gen);
npv = length(pv);
npq = length(pq);
npvq = npv + npq;

%% construct initial system v1, v2, zr, zi, v, z, x
t = 'construct initial system v, z';
sv1 = ac.params_var('va');
sv2 = ac.params_var('vm');
szr = ac.params_var('zr');
szi = ac.params_var('zi');

%% do one Newton step to get voltages that are not solution but not flat start
mpopt1 = mpoption(mpopt, 'pf.nr.max_it', 1);
[xpf, success, i] = ac.solve_power_flow(mpc, mpopt1);
sv1([pv; pq]) = xpf(1:npvq);
sv2(pq)       = xpf(npvq+(1:npq));

%% adjust values of z
szr(1) = 0.67;
szi = [0.1; 0.2; 0.3];

%% initialize v, z, x
sv = sv2 .* exp(1j * sv1);
sz = szr + 1j * szi;
sx = [sv; sz];
nx = length(sx);
t_is(nx, nv+nz, 12, t);

%%-----  tests using system voltages  -----
t = 'ac.port_inj_power(x) : ';
v10 = sv1; v20 = sv2; zr0 = szr; zi0 = szi; %% init w/ system v, z components
v0 = v20 .* exp(1j * v10);
z0 = zr0 + 1j * zi0;
x0 = [v0; z0];
Nv = length(v0);
Nz = length(z0);
[S0, Sv1, Sv2, Szr, Szi] = ac.port_inj_power(x0);   %% analytical

%% check matrix input/output
SS0 = ac.port_inj_power(x0*ones(1,Nv));
t_is(SS0, S0*ones(1,Nv), 12, [t 'matrix input']);

%% Sv1
v = (v20*ones(1,Nv)) .* exp(1j * (v10*ones(1,Nv) + dx*eye(Nv,Nv)));
z = (zr0 + 1j * zi0) * ones(1,Nv);
x = [v; z];
SS = ac.port_inj_power(x);
num_Sv1b = (SS - SS0) / dx;
t_is(full(Sv1), num_Sv1b, 6, [t 'Sv1']);

%% Sv2
v = (v20*ones(1,Nv) + dx*eye(Nv,Nv)) .* exp(1j * (v10*ones(1,Nv)));
z = (zr0 + 1j * zi0) * ones(1,Nv);
x = [v; z];
SS = ac.port_inj_power(x);
num_Sv2b = (SS - SS0) / dx;
t_is(full(Sv2), num_Sv2b, 6, [t 'Sv2']);

SS0 = ac.port_inj_power(x0*ones(1,Nz));

%% Szr
v = (v20*ones(1,Nz)) .* exp(1j * (v10*ones(1,Nz)));
z = (zr0 * ones(1,Nz) + dx*eye(Nz,Nz)) + 1j * (zi0 * ones(1,Nz));
x = [v; z];
SS = ac.port_inj_power(x);
num_Szrb = (SS - SS0) / dx;
t_is(full(Szr), num_Szrb, 6, [t 'Szr']);

%% Szi
v = (v20*ones(1,Nz)) .* exp(1j * (v10*ones(1,Nz)));
z = (zr0 * ones(1,Nz)) + 1j * (zi0 * ones(1,Nz) + dx*eye(Nz,Nz));
x = [v; z];
SS = ac.port_inj_power(x);
num_Szib = (SS - SS0) / dx;
t_is(full(Szi), num_Szib, 6, [t 'Szi']);

t = 'ac.port_inj_power(x, 1, idx) : ';
[iS0, iSv1, iSv2, iSzr, iSzi] = ac.port_inj_power(x0, 1, idx);
t_is(iS0, S0(idx), 12, [t 'S0']);
t_is(iSv1, Sv1(idx, :), 12, [t 'Sv1']);
t_is(iSv2, Sv2(idx, :), 12, [t 'Sv2']);
t_is(iSzr, Szr(idx, :), 12, [t 'Szr']);
t_is(iSzi, Szi(idx, :), 12, [t 'Szi']);

t = 'ac.port_inj_power_hess(x0, ek) == ac.p_i_p_h(x0, 1, 1, k) : ';
for k = 1:length(lam)
    ek = e0; ek(k) = 1;
    H1 = ac.port_inj_power_hess(x0, ek);
    H2 = ac.port_inj_power_hess(x0, 1, 1, k);
    t_is(H1, H2, 12, sprintf('%s%d', t, k));
end

t = 'ac.port_inj_power_hess(x, lam) : ';
H = ac.port_inj_power_hess(x0, lam);
HH = sparse(size(H, 1), size(H, 2));
for k = 1:length(lam)
    ek = e0; ek(k) = 1;
    HH = HH + lam(k) * ac.port_inj_power_hess(x0, ek);
end
t_is(H, HH, 12, [t 'weighted sum of indiv Hessians']);

t = 'ac.port_inj_power_hess(x, lam, 1, idx) : ';
H = ac.port_inj_power_hess(x0, lam(idx), 1, idx);
HH = sparse(size(H, 1), size(H, 2));
for k = 1:length(idx)
    ek = e0; ek(idx(k)) = 1;
    HH = HH + lam(idx(k)) * ac.port_inj_power_hess(x0, ek);
end
t_is(H, HH, 12, [t 'weighted sum of indiv Hessians']);

t = 'ac.port_inj_power_hess(x, lam) : ';
H = ac.port_inj_power_hess(x0, lam);
[S0, Sv1, Sv2, Szr, Szi] = ac.port_inj_power(x0);
numH = zeros(2*nx, 2*nx);
for k = 1:Nv
    v1 = v10; v1(k) = v1(k) + dx;
    v = v20 .* exp(1j * v1);
    x = [v; z0];
    [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x);
    numH(:, k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;

    v2 = v20; v2(k) = v2(k) + dx;
    v = v2 .* exp(1j * v10);
    x = [v; z0];
    [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x);
    numH(:, Nv+k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;
end
for k = 1:Nz
    z = zr0 + 1j * zi0;
    z(k) = z(k) + dx;
    x = [v0; z];
    [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x);
    numH(:, 2*Nv+k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;

    z = zr0 + 1j * zi0;
    z(k) = z(k) + 1j*dx;
    x = [v0; z];
    [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x);
    numH(:, 2*Nv+Nz+k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;
end
t_is(full(H), numH, 5, [t 'numerical Hessian']);

%%-----  tests using port voltages  -----
t = 'ac.port_inj_power(x, 0) : ';
v10 = C'*sv1; v20 = C'*sv2; zr0 = D'*szr; zi0 = D'*szi; %% init w/ port v, z components
v0 = v20 .* exp(1j * v10);
z0 = zr0 + 1j * zi0;
x0 = [v0; z0];
Nv = length(v0);
Nz = length(z0);
[S0, Sv1, Sv2, Szr, Szi] = ac.port_inj_power(x0, 0);    %% analytical

%% check matrix input/output
SS0 = ac.port_inj_power(x0*ones(1,Nv), 0);
t_is(SS0, S0*ones(1,Nv), 12, [t 'matrix input']);

%% Sv1
v = (v20*ones(1,Nv)) .* exp(1j * (v10*ones(1,Nv) + dx*eye(Nv,Nv)));
z = (zr0 + 1j * zi0) * ones(1,Nv);
x = [v; z];
SS = ac.port_inj_power(x, 0);
num_Sv1b = (SS - SS0) / dx;
t_is(full(Sv1), num_Sv1b, 6, [t 'Sv1']);

%% Sv2
v = (v20*ones(1,Nv) + dx*eye(Nv,Nv)) .* exp(1j * (v10*ones(1,Nv)));
z = (zr0 + 1j * zi0) * ones(1,Nv);
x = [v; z];
SS = ac.port_inj_power(x, 0);
num_Sv2b = (SS - SS0) / dx;
t_is(full(Sv2), num_Sv2b, 6, [t 'Sv2']);

SS0 = ac.port_inj_power(x0*ones(1,Nz), 0);

%% Szr
v = (v20*ones(1,Nz)) .* exp(1j * (v10*ones(1,Nz)));
z = (zr0 * ones(1,Nz) + dx*eye(Nz,Nz)) + 1j * (zi0 * ones(1,Nz));
x = [v; z];
SS = ac.port_inj_power(x, 0);
num_Szrb = (SS - SS0) / dx;
t_is(full(Szr), num_Szrb, 6, [t 'Szr']);

%% Szi
v = (v20*ones(1,Nz)) .* exp(1j * (v10*ones(1,Nz)));
z = (zr0 * ones(1,Nz)) + 1j * (zi0 * ones(1,Nz) + dx*eye(Nz,Nz));
x = [v; z];
SS = ac.port_inj_power(x, 0);
num_Szib = (SS - SS0) / dx;
t_is(full(Szi), num_Szib, 6, [t 'Szi']);

t = 'ac.port_inj_power(x, 0, idx) : ';
[iS0, iSv1, iSv2, iSzr, iSzi] = ac.port_inj_power(x0, 0, idx);
t_is(iS0, S0(idx), 12, [t 'S0']);
t_is(iSv1, Sv1(idx, :), 12, [t 'Sv1']);
t_is(iSv2, Sv2(idx, :), 12, [t 'Sv2']);
t_is(iSzr, Szr(idx, :), 12, [t 'Szr']);
t_is(iSzi, Szi(idx, :), 12, [t 'Szi']);

t = 'ac.port_inj_power_hess(x0, ek, 0) == ac.p_i_p_h(x0, 1, 0, k) : ';
for k = 1:length(lam)
    ek = e0; ek(k) = 1;
    H1 = ac.port_inj_power_hess(x0, ek, 0);
    H2 = ac.port_inj_power_hess(x0, 1, 0, k);
    t_is(H1, H2, 12, sprintf('%s%d', t, k));
end

t = 'ac.port_inj_power_hess(x, lam, 0) : ';
H = ac.port_inj_power_hess(x0, lam, 0);
HH = sparse(size(H, 1), size(H, 2));
for k = 1:length(lam)
    ek = e0; ek(k) = 1;
    HH = HH + lam(k) * ac.port_inj_power_hess(x0, ek, 0);
end
t_is(H, HH, 12, [t 'weighted sum of indiv Hessians']);

t = 'ac.port_inj_power_hess(x, lam, 0, idx) : ';
H = ac.port_inj_power_hess(x0, lam(idx), 0, idx);
HH = sparse(size(H, 1), size(H, 2));
for k = 1:length(idx)
    ek = e0; ek(idx(k)) = 1;
    HH = HH + lam(idx(k)) * ac.port_inj_power_hess(x0, ek, 0);
end
t_is(H, HH, 12, [t 'weighted sum of indiv Hessians']);

t = 'ac.port_inj_power_hess(x, lam, 0) : ';
H = ac.port_inj_power_hess(x0, lam, 0);
[S0, Sv1, Sv2, Szr, Szi] = ac.port_inj_power(x0, 0);
Nx = 2*Nv+2*Nz;
numH = zeros(Nx, Nx);
for k = 1:Nv
    v1 = v10; v1(k) = v1(k) + dx;
    v = v20 .* exp(1j * v1);
    x = [v; z0];
    [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x, 0);
    numH(:, k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;

    v2 = v20; v2(k) = v2(k) + dx;
    v = v2 .* exp(1j * v10);
    x = [v; z0];
    [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x, 0);
    numH(:, Nv+k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;
end
for k = 1:Nz
    z = zr0 + 1j * zi0;
    z(k) = z(k) + dx;
    x = [v0; z];
    [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x, 0);
    numH(:, 2*Nv+k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;

    z = zr0 + 1j * zi0;
    z(k) = z(k) + 1j*dx;
    x = [v0; z];
    [S0p, Sv1p, Sv2p, Szrp, Szip] = ac.port_inj_power(x, 0);
    numH(:, 2*Nv+Nz+k) = ([Sv1p, Sv2p, Szrp, Szip]- [Sv1, Sv2, Szr, Szi]).' * lam / dx;
end
t_is(full(H), numH, 5, [t 'numerical Hessian']);


%  x0 = [va([pv; pq]); vm(pq)];
% pf_fcn = @(x)power_flow_equations(ac, x, va, vm, z, ref, pv, pq);
% [F, J] = power_flow_equations(obj, x, va, vm, z, ref, pv, pq)
% [F, J] = pf_fcn(ac, x)




% %% analytical
% [g0, dg] = ac.opf_power_balance_fcn(x0);
% 
% %% compute Jacobian numerically via finite differences
% num_dg = zeros(size(dg));
% nx = length(x0);
% dx = 1e-8;
% for j = 1:nx
%     x = x0;
%     x(j) = x(j) + dx;
%     g = ac.opf_power_balance_fcn(x);
%     num_dg(:, j) = (g - g0) / dx;
% end
% 
% t = '[g, dg] = ac.opf_power_balance_fcn(x) : dg';
% t_is(full(dg), num_dg, 6, t);
% 
% full(dg)
% num_dg


% disp(ac)

t_end;

if nargout
    obj = ac;
end

return;
