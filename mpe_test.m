define_constants;
mpopt = mpoption('out.all', 0, 'verbose', 2);

% mpc = rundcpf(loadcase('t_case9_opfv2'), mpopt);
% dc = mpe_network_dc().create_model(mpc);
% 
% mpc = ext2int(runpf(loadcase('t_case9_opfv2'), mpopt));
% ac = mpe_network_acps().create_model(mpc);

% mpc = ext2int(loadcase('t_case9_opfv2'));
mpc = ext2int(loadcase('t_case9_gizmo'));
% mpc = ext2int(loadcase('case9'));
% mpc = rmfield(mpc, 'order');

% mpc = loadcase('case9');
% mpc.branch(4, RATE_A) = 80;
% mpc = ext2int(mpc);

ac = mpe_network_acps().create_model(mpc);

[x, success, i] = ac.solve_power_flow(mpc, mpopt);
[x, success, i] = ac.solve_opf(mpc, mpopt);


% mpc = ext2int(loadcase('t_case9_gizmo'));
% mpc = rmfield(mpc, 'order');
ac = mpe_network_acps_test().create_model(mpc);

[x, success, i] = ac.solve_power_flow(mpc, mpopt);
[x, success, i] = ac.solve_opf(mpc, mpopt);


mpc = ext2int(loadcase('t_case9_opf'));

dc = mpe_network_dc().create_model(mpc);
[x, success, i] = dc.solve_opf(mpc, mpopt);
