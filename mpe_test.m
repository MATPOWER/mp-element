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
norm(P-P1)
dc.C{1} * P
gen = dc.mpe_by_name('gen');
Pg = gen.port_inj_power(x, 1)
Pg31 = gen.port_inj_power(x, 1, [3;1])
Pg2 = gen.port_inj_power(x, 1, 2)

% return;

mpc = ext2int(runpf(loadcase('t_case9_opfv2'), mpopt));
mpc = ext2int(loadcase('t_case9_opfv2'));
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
ac.C{1} * I
gen = ac.mpe_by_name('gen');
Sg = gen.port_inj_power(x, 1)
Sg31 = gen.port_inj_power(x, 1, [3;1])
Sg2 = gen.port_inj_power(x, 1, 2)

[S, dSdva, dSdvm, dSdzr, dSdzi] = ac.port_inj_power(x, 1, [3;2;1])

[x, success, i] = ac.solve_power_flow(mpc)

% [S, dSdva, dSdvm, dSdzr, dSdzi] = gen.port_inj_power(x, 1)
% C = horzcat(gen.C{:});
% D = horzcat(gen.D{:});
% A = [   C sparse(size(C, 1), size(D, 2));
%         sparse(size(D, 1), size(C, 2)) D   ];
% [S, dSdva, dSdvm, dSdzr, dSdzi] = gen.port_inj_power(A'*x, 0)


% mpe_type_lib = struct( ...
%     'acsp', {{ @acsp_gen, @acsp_load, @acsp_branch }}, ...
%     'dc',   {{ @dc_gen, @dc_load, @dc_branch }} ...
% );
% % mpe_types = { @mp_gen, @mp_load, @mp_branch };
% formulation = 'acsp';
% % formulation = 'dc';
% mpe_types = mpe_type_lib.(formulation);
% 
% %% bus to node mappings
% nn = size(mpc.bus, 1);
% n2b = mpc.bus(:, BUS_I);
% b2n = sparse(max(n2b), 1);
% b2n(n2b) = (1:nn)';
% nb = nn;
% 
% mpc0 = mpc;
% mpc = ext2int(mpc0);
% 
% om = opt_model();
% 
% %% set up bus voltage state variables
% refs = find(mpc.bus(:, BUS_TYPE) == REF);
% Vm   = mpc.bus(:, VM);
% Va   = mpc.bus(:, VA) * (pi/180);
% Vau = Inf(nb, 1);       %% voltage angle limits
% Val = -Vau;
% Vau(refs) = Va(refs);   %% voltage angle reference constraints
% Val(refs) = Va(refs);
% switch formulation
%     case 'acsp'
%         om.add_var('Va', nb, Va, Val, Vau);
%         om.add_var('Vm', nb, Vm, mpc.bus(:, VMIN), mpc.bus(:, VMAX));
%     case 'dc'
%         om.add_var('Va', nb, Va, Val, Vau);
%     otherwise
%         error('formulation ''%s'' not implemented', formulation)
% end
% 
% %% create relevant element objects
% mpe_list = {};
% for mpe_type = mpe_types
%     %% construct element type
%     mpe = mpe_type{1}(formulation, mpc);
%     if mpe.count(mpc)
%         mpe_list{end+1} = mpe;
%     end
% end
% 
% %% add vars
% for mpe = mpe_list
%     om = mpe{1}.add_vars(om, mpc);
% end
% 
% %% build params
% for mpe = mpe_list
%     mpe{1}.build_params(om, mpc);
% end
% 
% 
% 
% %% display
% mpe_list{:}
% om
% for t = 1:length(mpe_list)
%     mpe = mpe_list{t};
%     mpe.name
%     fields = mpe.model_params();
%     for f = 1:length(fields)
%         if ~isempty(mpe.(fields{f}))
%             fields{f}
%             mpe.(fields{f})
%         end
%     end
% end
