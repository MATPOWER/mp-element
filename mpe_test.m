define_constants;
mpopt = mpoption('out.all', 0, 'verbose', 2);
mpc = rundcpf(loadcase('t_case9_opfv2'), mpopt);
nb = size(mpc.bus, 1);
ng = size(mpc.gen, 1);

dc = dc_aggregate();
dc.create_model(mpc);
nv = dc.nv;
nz = dc.nz;
np = dc.np;

dc
v = mpc.bus(:, VA) * pi/180;
z = mpc.gen(:, PG) / mpc.baseMVA;
x = [v;z];
A = [   dc.C{1} sparse(nv, nz);
        sparse(nz, np) dc.D{1}   ];
P  = dc.port_inj_power(x, 1)
P1 = dc.port_inj_power(A'*x, 0)
P2  = dc.port_inj_power(x, 1, [3;1])
norm(P-P1)
norm(P([3;1])-P2)
dc.C{1} * P
gen = dc.mpe_by_name('gen');
Pg = gen.port_inj_power(x, 1)
Pg31 = gen.port_inj_power(x, 1, [3;1])
Pg2 = gen.port_inj_power(x, 1, 2)

% return;

mpc = ext2int(runpf(loadcase('t_case9_opfv2'), mpopt));
% mpc = ext2int(loadcase('t_case9_opfv2'));
mpc = rmfield(mpc, 'order');
ac = acsp_aggregate();
ac.create_model(mpc);

ac
v = mpc.bus(:, VM) .* exp(1j * mpc.bus(:, VA) * pi/180);
z = (mpc.gen(:, PG) + 1j * mpc.gen(:, QG)) / mpc.baseMVA;
x = [v;z];
A = [   dc.C{1} sparse(nv, nz);
        sparse(nz, np) dc.D{1}   ];
S  = ac.port_inj_power(x, 1)
S1 = ac.port_inj_power(A'*x, 0)
norm(S-S1)
ac.C{1} * S
I  = ac.port_inj_current(x, 1)
I1 = ac.port_inj_current(A'*x, 0)
norm(I-I1)
ac.C{1} * I

[I, Iva, Ivm, Izr, Izi] = ac.port_inj_current(x, 1, [3;2;1])
[I, Iva, Ivm, Izr, Izi] = ac.port_inj_current(x, 1)

% return;

gen = ac.mpe_by_name('gen');
Sg = gen.port_inj_power(x, 1)
Sg31 = gen.port_inj_power(x, 1, [3;1])
Sg2 = gen.port_inj_power(x, 1, 2)

[S, Sva, Svm, Szr, Szi] = ac.port_inj_power(x, 1, [3;2;1])
[S, Sva, Svm, Szr, Szi] = ac.port_inj_power(x, 1)

np = ac.np;
lam = [1:ac.np]' / 100;

H = ac.port_inj_power_hess(x, lam, 1)
H = ac.port_inj_power_hess(x, lam(1:3), 1, [3;2;1])

mpc = ext2int(loadcase('t_case9_opfv2'));
mpc = rmfield(mpc, 'order');
ac = acsp_aggregate();
ac.create_model(mpc);
[x, success, i] = ac.solve_power_flow(mpc)

[x, success, i] = ac.solve_opf(mpc)
