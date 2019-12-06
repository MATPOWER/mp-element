define_constants;
mpopt = mpoption('out.all', 0, 'verbose', 2);

% mpc = rundcpf(loadcase('t_case9_opfv2'), mpopt);
% dc = dc_aggregate().create_model(mpc);
% 
% mpc = ext2int(runpf(loadcase('t_case9_opfv2'), mpopt));
% ac = acsp_aggregate().create_model(mpc);

mpc = ext2int(loadcase('t_case9_opfv2'));
mpc = rmfield(mpc, 'order');
ac = acsp_aggregate().create_model(mpc)

[x, success, i] = ac.solve_power_flow(mpc, mpopt)

[x, success, i] = ac.solve_opf(mpc)
