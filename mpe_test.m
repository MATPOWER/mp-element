define_constants;
mpc = loadcase('t_case9_opfv2');

dc = dc_aggregate();
dc.create_model(mpc);

ac = acsp_aggregate();
ac.create_model(mpc);

dc
dc.node
dc.var
{dc.var.order.name}
dc.mpe_list{1}
dc.mpe_list{2}
dc.mpe_list{3}
dc.mpe_list{4}

ac
ac.node
ac.var
{ac.var.order.name}
ac.mpe_list{1}
ac.mpe_list{2}
ac.mpe_list{3}
ac.mpe_list{4}



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
