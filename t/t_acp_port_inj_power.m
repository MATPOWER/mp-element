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

t_begin(10, quiet);

define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

casefile = 't_case9_opfv2';

mpopt = mpoption('out.all', 0, 'verbose', 0);

dx = 1e-8;

%% create aggregate system object
mpc = ext2int(loadcase('t_case9_opfv2'));
mpc = rmfield(mpc, 'order');
ac = acsp_aggregate().create_model(mpc);
C = ac.getC();
D = ac.getD();
np = ac.np;
nv = ac.nv/2;
nz = ac.nz;
A = [   C sparse(nv, nz);
        sparse(nz, np) D    ];

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

t = 'ac.port_inj_power(x) : ';
v10 = sv1; v20 = sv2; zr0 = szr; zi0 = szi; %% init w/ system v, z components
v0 = v20 .* exp(1j * v10);
z0 = zr0 + 1j * zi0;
x0 = [v0; z0];
Nv = length(v0);
Nz = length(z0);
[S0, Sv1, Sv2, Szr, Szi] = ac.port_inj_power(x0);   %% analytical

% %% compute Jacobian numerically via finite differences
% num_Sv1 = zeros(size(Sv1));
% num_Sv2 = zeros(size(Sv2));
% num_Szr = zeros(size(Szr));
% num_Szi = zeros(size(Szi));
% for j = 1:Nv
%     v1 = v10; v2 = v20; zr = zr0; zi = zi0;
%     v1(j) = v1(j) + dx;     %% perturbation
%     v = v2 .* exp(1j * v1); z = zr + 1j * zi; x = [v; z];
%     S = ac.port_inj_power(x);
%     num_Sv1(:, j) = (S - S0) / dx;
% end
% t_is(full(Sv1), num_Sv1, 6, [t 'Sv1']);
% 
% for j = 1:Nv
%     v1 = v10; v2 = v20; zr = zr0; zi = zi0;
%     v2(j) = v2(j) + dx;     %% perturbation
%     v = v2 .* exp(1j * v1); z = zr + 1j * zi; x = [v; z];
%     S = ac.port_inj_power(x);
%     num_Sv2(:, j) = (S - S0) / dx;
% end
% t_is(full(Sv2), num_Sv2, 6, [t 'Sv2']);
% 
% for j = 1:Nz
%     v1 = v10; v2 = v20; zr = zr0; zi = zi0;
%     zr(j) = zr(j) + dx;     %% perturbation
%     v = v2 .* exp(1j * v1); z = zr + 1j * zi; x = [v; z];
%     S = ac.port_inj_power(x);
%     num_Szr(:, j) = (S - S0) / dx;
% end
% t_is(full(Szr), num_Szr, 6, [t 'Szr']);
% 
% for j = 1:Nz
%     v1 = v10; v2 = v20; zr = zr0; zi = zi0;
%     zi(j) = zi(j) + dx;     %% perturbation
%     v = v2 .* exp(1j * v1); z = zr + 1j * zi; x = [v; z];
%     S = ac.port_inj_power(x);
%     num_Szi(:, j) = (S - S0) / dx;
% end
% t_is(full(Szi), num_Szi, 6, [t 'Szi']);

%% check matrix input/output
SS0 = ac.port_inj_power(x0*ones(1,Nv));
t_is(SS0, S0*ones(1,Nv), 12, [t 'matrix input']);

t = 'ac.port_inj_power(x) : ';
v = (v20*ones(1,Nv)) .* exp(1j * (v10*ones(1,Nv) + dx*eye(Nv,Nv)));
z = (zr0 + 1j * zi0) * ones(1,Nv);
x = [v; z];
SS = ac.port_inj_power(x);
num_Sv1b = (SS - SS0) / dx;
t_is(full(Sv1), num_Sv1b, 6, [t 'Sv1']);

v = (v20*ones(1,Nv) + dx*eye(Nv,Nv)) .* exp(1j * (v10*ones(1,Nv)));
z = (zr0 + 1j * zi0) * ones(1,Nv);
x = [v; z];
SS = ac.port_inj_power(x);
num_Sv2b = (SS - SS0) / dx;
t_is(full(Sv2), num_Sv2b, 6, [t 'Sv2']);

SS0 = ac.port_inj_power(x0*ones(1,Nz));

v = (v20*ones(1,Nz)) .* exp(1j * (v10*ones(1,Nz)));
z = (zr0 * ones(1,Nz) + dx*eye(Nz,Nz)) + 1j * (zi0 * ones(1,Nz));
x = [v; z];
SS = ac.port_inj_power(x);
num_Szrb = (SS - SS0) / dx;
t_is(full(Szr), num_Szrb, 6, [t 'Szr']);

v = (v20*ones(1,Nz)) .* exp(1j * (v10*ones(1,Nz)));
z = (zr0 * ones(1,Nz)) + 1j * (zi0 * ones(1,Nz) + dx*eye(Nz,Nz));
x = [v; z];
SS = ac.port_inj_power(x);
num_Szib = (SS - SS0) / dx;
t_is(full(Szi), num_Szib, 6, [t 'Szi']);

t = 'ac.port_inj_power(x, 0) : ';
v10 = C'*sv1; v20 = C'*sv2; zr0 = D'*szr; zi0 = D'*szi; %% init w/ port v, z components
v0 = v20 .* exp(1j * v10);
z0 = zr0 + 1j * zi0;
x0 = [v0; z0];
Nv = length(v0);
Nz = length(z0);
[S0, Sv1, Sv2, Szr, Szi] = ac.port_inj_power(x0, 0);    %% analytical

%% compute Jacobian numerically via finite differences
num_Sv1 = zeros(size(Sv1));
num_Sv2 = zeros(size(Sv2));
num_Szr = zeros(size(Szr));
num_Szi = zeros(size(Szi));
dx = 1e-8;
for j = 1:Nv
    v1 = v10; v2 = v20; zr = zr0; zi = zi0;
    v1(j) = v1(j) + dx;     %% perturbation
    v = v2 .* exp(1j * v1); z = zr + 1j * zi; x = [v; z];
    S = ac.port_inj_power(x, 0);
    num_Sv1(:, j) = (S - S0) / dx;
end
t_is(full(Sv1), num_Sv1, 6, [t 'Sv1']);

for j = 1:Nv
    v1 = v10; v2 = v20; zr = zr0; zi = zi0;
    v2(j) = v2(j) + dx;     %% perturbation
    v = v2 .* exp(1j * v1); z = zr + 1j * zi; x = [v; z];
    S = ac.port_inj_power(x, 0);
    num_Sv2(:, j) = (S - S0) / dx;
end
t_is(full(Sv2), num_Sv2, 6, [t 'Sv2']);

for j = 1:Nz
    v1 = v10; v2 = v20; zr = zr0; zi = zi0;
    zr(j) = zr(j) + dx;     %% perturbation
    v = v2 .* exp(1j * v1); z = zr + 1j * zi; x = [v; z];
    S = ac.port_inj_power(x, 0);
    num_Szr(:, j) = (S - S0) / dx;
end
t_is(full(Szr), num_Szr, 6, [t 'Szr']);

for j = 1:Nz
    v1 = v10; v2 = v20; zr = zr0; zi = zi0;
    zi(j) = zi(j) + dx;     %% perturbation
    v = v2 .* exp(1j * v1); z = zr + 1j * zi; x = [v; z];
    S = ac.port_inj_power(x, 0);
    num_Szi(:, j) = (S - S0) / dx;
end
t_is(full(Szi), num_Szi, 6, [t 'Szi']);


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
