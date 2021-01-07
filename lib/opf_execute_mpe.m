function [results, success, raw] = opf_execute_mpe(opf, mpopt)

om = opf.mm;

%%-----  setup  -----
%% options
dc  = strcmp(upper(mpopt.model), 'DC');
alg = upper(mpopt.opf.ac.solver);
vcart = ~dc && mpopt.opf.v_cartesian;

%% get indexing
[vv, ll, nne, nni] = om.get_idx();

if mpopt.verbose > 0
    v = mpver('all');
    fprintf('\nMATPOWER Version %s, %s', v.Version, v.Date);
end

if dc
    %%-----  run DC OPF solver  -----
    if mpopt.verbose > 0
        fprintf(' -- DC Optimal Power Flow\n');
    end
    [results, success, raw] = dcopf_solver_mpe(opf, mpopt);
else
    %%-----  run AC OPF solver  -----
    if mpopt.verbose > 0
        fprintf(' -- AC Optimal Power Flow\n  AC OPF formulation: ');
        if vcart
            v = 'cartesian';
        else
            v = 'polar';
        end
        if mpopt.opf.current_balance
            v2 = 'current';
        else
            v2 = 'power';
        end
        fprintf('%s voltages, %s balance eqns\n', v, v2);
    end

    %% run specific AC OPF solver
    switch alg
        case 'IPOPT'
            if ~have_feature('ipopt')
                error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires IPOPT (see https://github.com/coin-or/Ipopt)', alg);
            end
        case 'FMINCON'
            if ~have_feature('fmincon')
                error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires FMINCON (Optimization Toolbox 2.x or later)', alg);
            end
        case 'KNITRO'
            if ~have_feature('knitro')
                error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires Artelys Knitro (see https://www.artelys.com/solvers/knitro/)', alg);
            end
    end
    [results, success, raw] = nlpopf_solver_mpe(opf, mpopt);
end %% if dc
if ~isfield(raw, 'output') || ~isfield(raw.output, 'alg') || isempty(raw.output.alg)
    raw.output.alg = alg;
end

if success
  if ~dc

    %% compute g, dg, f, df, d2f if requested by opf.return_raw_der = 1
    if mpopt.opf.return_raw_der
      %% move from results to raw if using v4.0 of MINOPF or TSPOPF
      if isfield(results, 'dg')
        raw.dg = results.dg;
        raw.g = results.g;
      end
      %% compute g, dg, unless already done by post-v4.0 MINOPF or TSPOPF
      if ~isfield(raw, 'dg')
        mpc = om.get_mpc();
        [g, geq, dg, dgeq] = nlp_consfcn(om, results.x);
        raw.g = [ geq; g];
        raw.dg = [ dgeq'; dg'];   %% true Jacobian organization
      end
      %% compute df, d2f
      [f, df, d2f] = nlp_costfcn(om, results.x);
      raw.df = df;
      raw.d2f = d2f;
    end
  end

  %% delete g and dg fields from results if using v4.0 of MINOPF or TSPOPF
  if isfield(results, 'dg')
    rmfield(results, 'dg');
    rmfield(results, 'g');
  end
else
  %% assign empty g, dg, f, df, d2f if requested by opf.return_raw_der = 1
  if ~dc && mpopt.opf.return_raw_der
    raw.dg = [];
    raw.g = [];
    raw.df = [];
    raw.d2f = [];
  end
end

%% assign values and limit shadow prices for variables
om_var_order = om.get('var', 'order');
for k = 1:length(om_var_order)
  name = om_var_order(k).name;
  if om.getN('var', name)
    idx = vv.i1.(name):vv.iN.(name);
    results.var.val.(name) = results.x(idx);
    results.var.mu.l.(name) = results.mu.var.l(idx);
    results.var.mu.u.(name) = results.mu.var.u(idx);
  end
end

%% assign shadow prices for linear constraints
om_lin_order = om.get('lin', 'order');
for k = 1:length(om_lin_order)
  name = om_lin_order(k).name;
  if om.getN('lin', name)
    idx = ll.i1.(name):ll.iN.(name);
    results.lin.mu.l.(name) = results.mu.lin.l(idx);
    results.lin.mu.u.(name) = results.mu.lin.u(idx);
  end
end

%% assign shadow prices for nonlinear constraints
if ~dc
  om_nle_order = om.get('nle', 'order');
  for k = 1:length(om_nle_order)
    name = om_nle_order(k).name;
    if om.getN('nle', name)
      results.nle.lambda.(name) = results.mu.nle(nne.i1.(name):nne.iN.(name));
    end
  end

  om_nli_order = om.get('nli', 'order');
  for k = 1:length(om_nli_order)
    name = om_nli_order(k).name;
    if om.getN('nli', name)
      results.nli.mu.(name) = results.mu.nli(nni.i1.(name):nni.iN.(name));
    end
  end
end

%% assign values for components of quadratic cost
om_qdc_order = om.get('qdc', 'order');
for k = 1:length(om_qdc_order)
  name = om_qdc_order(k).name;
  if om.getN('qdc', name)
    results.qdc.(name) = om.eval_quad_cost(results.x, name);
  end
end

%% assign values for components of general nonlinear cost
om_nlc_order = om.get('nlc', 'order');
for k = 1:length(om_nlc_order)
  name = om_nlc_order(k).name;
  if om.getN('nlc', name)
    results.nlc.(name) = om.eval_nln_cost(results.x, name);
  end
end

%% assign values for components of legacy user cost
om_cost_order = om.get('cost', 'order');
for k = 1:length(om_cost_order)
  name = om_cost_order(k).name;
  if om.getN('cost', name)
    results.cost.(name) = om.eval_legacy_cost(results.x, name);
  end
end

%% if single-block PWL costs were converted to POLY, insert dummy y into x
%% Note: The "y" portion of x will be nonsense, but everything should at
%%       least be in the expected locations.
pwl1 = om.get_userdata('pwl1');
if ~isempty(pwl1)
  %% get indexing
  vv = om.get_idx();
  if dc
    nx = vv.iN.Pg;
  else
    nx = vv.iN.Qg;
  end
  y = zeros(length(pwl1), 1);
  raw.xr = [ raw.xr(1:nx); y; raw.xr(nx+1:end)];
  results.x = [ results.x(1:nx); y; results.x(nx+1:end)];
end
