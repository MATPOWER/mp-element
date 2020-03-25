define_constants;
mpopt = mpoption('out.all', 0, 'verbose', 2);

% mpc = rundcpf(loadcase('t_case9_opfv2'), mpopt);
% dc = dc_aggregate().create_model(mpc);
% 
% mpc = ext2int(runpf(loadcase('t_case9_opfv2'), mpopt));
% ac = acsp_aggregate().create_model(mpc);

% mpc = ext2int(loadcase('t_case9_opfv2'));
mpc = ext2int(loadcase('t_case9_gizmo'));
% mpc = ext2int(loadcase('case9'));
% mpc = rmfield(mpc, 'order');

% mpc = loadcase('case9');
% mpc.branch(4, RATE_A) = 80;
% mpc = ext2int(mpc);

ac = acsp_aggregate().create_model(mpc);

[x, success, i] = ac.solve_power_flow(mpc, mpopt);
[x, success, i] = ac.solve_opf(mpc, mpopt);


% mpc = ext2int(loadcase('t_case9_gizmo'));
% mpc = rmfield(mpc, 'order');
ac = acsp_test_aggregate().create_model(mpc);

[x, success, i] = ac.solve_power_flow(mpc, mpopt);
[x, success, i] = ac.solve_opf(mpc, mpopt);


mpc = ext2int(loadcase('t_case9_opf'));

dc = dc_aggregate().create_model(mpc);
[x, success, i] = dc.solve_opf(mpc, mpopt);
