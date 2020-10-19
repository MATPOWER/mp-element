function om = opf_setup_mpe(mpc, mpopt)
%OPF  Constructs an OPF model object from a MATPOWER case struct.
%   OM = OPF_SETUP(MPC, MPOPT)
%
%   Assumes that MPC is a MATPOWER case struct with internal indexing,
%   all equipment in-service, etc.
%
%   See also OPF, EXT2INT, OPF_EXECUTE.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if mpopt.verbose, fprintf('-----  opf_setup() : MPE version  -----\n'); end

%% options
dc  = strcmp(upper(mpopt.model), 'DC');
alg = upper(mpopt.opf.ac.solver);
use_vg = mpopt.opf.use_vg;
vcart = ~dc && mpopt.opf.v_cartesian;

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% define flag to indicate whether we are tied to legacy formulation
%% implemented by legacy MINOS and PDIPM solvers (e.g. with hard-coded
%% costs and constrain
if strcmp(alg, 'MINOPF') || strcmp(alg, 'PDIPM') || ...
        strcmp(alg, 'TRALM') || strcmp(alg, 'SDPOPF')
    legacy_formulation = 1;
    if vcart
        error('Option ''opf.v_cartesian'' = 1 is not compatible with ''opf.solver.ac''=''%s''.', alg);
    end
    if mpopt.opf.current_balance
        error('Option ''opf.current_balance'' = 1 is not compatible with ''opf.solver.ac''=''%s''.', alg);
    end
else
    legacy_formulation = 0;
end
if legacy_formulation && ~dc
    error('opf_setup_mpe not compatible with ''opf.solver.ac''=''%s''.', alg);
end
if ~dc && ( ~isempty(mpopt.exp.sys_wide_zip_loads.pw) && ...
                    ~isequal(mpopt.exp.sys_wide_zip_loads.pw, [1 0 0]) || ...
            ~isempty(mpopt.exp.sys_wide_zip_loads.qw) && ...
                    ~isequal(mpopt.exp.sys_wide_zip_loads.qw, [1 0 0]) )
    if vcart
        warning('Voltage dependent loads are not supported with option ''opf.v_cartesian'' = 1. Reverting to constant power load model.');
    end
    if mpopt.opf.current_balance
        warning('Voltage dependent loads are not supported with option ''opf.current_balance'' = 1. Reverting to constant power load model.');
    end
end

%% data dimensions
nb   = size(mpc.bus, 1);    %% number of buses
ng   = size(mpc.gen, 1);    %% number of dispatchable injections
if isfield(mpc, 'A')
  nlin = size(mpc.A, 1);    %% number of linear user constraints
else
  nlin = 0;
end
if isfield(mpc, 'N')
  nw = size(mpc.N, 1);      %% number of general cost vars, w
else
  nw = 0;
end

if dc
  %% ignore reactive costs for DC
  mpc.gencost = pqcost(mpc.gencost, ng);

  %% reduce A and/or N from AC dimensions to DC dimensions, if needed
  if nlin || nw
    acc = [nb+(1:nb) 2*nb+ng+(1:ng)];   %% Vm and Qg columns
    if nlin && size(mpc.A, 2) >= 2*nb + 2*ng
      %% make sure there aren't any constraints on Vm or Qg
      if any(any(mpc.A(:, acc)))
        error('opf_setup: attempting to solve DC OPF with user constraints on Vm or Qg');
      end
      mpc.A(:, acc) = [];               %% delete Vm and Qg columns
    end
    if nw && size(mpc.N, 2) >= 2*nb + 2*ng
      %% make sure there aren't any costs on Vm or Qg
      if any(any(mpc.N(:, acc)))
        [ii, jj] = find(mpc.N(:, acc));
        ii = unique(ii);    %% indices of w with potential non-zero cost terms from Vm or Qg
        if any(mpc.Cw(ii)) || (isfield(mpc, 'H') && ~isempty(mpc.H) && ...
                any(any(mpc.H(:, ii))))
          error('opf_setup: attempting to solve DC OPF with user costs on Vm or Qg');
        end
      end
      mpc.N(:, acc) = [];               %% delete Vm and Qg columns
    end
  end
else    %% AC
  if use_vg     %% adjust bus voltage limits based on generator Vg setpoint
    %% gen connection matrix, element i, j is 1 if, generator j at bus i is ON
    Cg = sparse(mpc.gen(:, GEN_BUS), (1:ng)', mpc.gen(:, GEN_STATUS) > 0, nb, ng);
    Vbg = Cg * sparse(1:ng, 1:ng, mpc.gen(:, VG), ng, ng);
    Vmax = max(Vbg, [], 2); %% zero for non-gen buses, else max Vg of gens @ bus
    ib = find(Vmax);                %% buses with online gens
    Vmin = max(2*Cg - Vbg, [], 2);  %% same as Vmax, except min Vg of gens @ bus
    Vmin(ib) = 2 - Vmin(ib);

    if use_vg == 1      %% use Vg setpoint directly
        mpc.bus(ib, VMAX) = Vmax(ib);   %% max set by max Vg @ bus
        mpc.bus(ib, VMIN) = Vmin(ib);   %% min set by min Vg @ bus
        mpc.bus(ib, VM) = mpc.bus(ib, VMAX);
    elseif use_vg > 0 && use_vg < 1     %% fractional value
        %% use weighted avg between original Vmin/Vmax limits and Vg
        mpc.bus(ib, VMAX) = (1-use_vg) * mpc.bus(ib, VMAX) + use_vg * Vmax(ib);
        mpc.bus(ib, VMIN) = (1-use_vg) * mpc.bus(ib, VMIN) + use_vg * Vmin(ib);
    else
        error('opf_setup: option ''opf.use_vg'' (= %g) cannot be negative or greater than 1', use_vg);
    end
  end
end

%% convert single-block piecewise-linear costs into linear polynomial cost
pwl1 = find(mpc.gencost(:, MODEL) == PW_LINEAR & mpc.gencost(:, NCOST) == 2);
% p1 = [];
if ~isempty(pwl1)
  x0 = mpc.gencost(pwl1, COST);
  y0 = mpc.gencost(pwl1, COST+1);
  x1 = mpc.gencost(pwl1, COST+2);
  y1 = mpc.gencost(pwl1, COST+3);
  m = (y1 - y0) ./ (x1 - x0);
  b = y0 - m .* x0;
  mpc.gencost(pwl1, MODEL) = POLYNOMIAL;
  mpc.gencost(pwl1, NCOST) = 2;
  mpc.gencost(pwl1, COST:COST+1) = [m b];
end

%% warn if there is more than one reference bus
refs = find(mpc.bus(:, BUS_TYPE) == REF);
if length(refs) > 1 && mpopt.verbose > 0
  errstr = ['\nopf_setup: Warning: Multiple reference buses.\n', ...
              '           For a system with islands, a reference bus in each island\n', ...
              '           may help convergence, but in a fully connected system such\n', ...
              '           a situation is probably not reasonable.\n\n' ];
  fprintf(errstr);
end


%% choose network model class
if dc
    netmodel = @mpe_network_dc;
else
    if vcart
        if mpopt.opf.current_balance
            netmodel = @mpe_network_acci;
        else
%             netmodel = @mpe_network_accs_test_nln;
%             netmodel = @mpe_network_accs_test;
            netmodel = @mpe_network_accs;
        end
    else
        if mpopt.opf.current_balance
            netmodel = @mpe_network_acpi;
        else
%             netmodel = @mpe_network_acps_test_nln;
%             netmodel = @mpe_network_acps_test;
            netmodel = @mpe_network_acps;
        end
    end
end

%%-----  build network model from data model (mpc)  -----
% nm = netmodel().create_model(mpc);
nm = netmodel().create_model(mpc, mpopt);   %% remove mpopt once we get all data (e.g. exp.sys_wide_zip_loads) out of options
om = nm.setup_opf(mpc, mpopt);

%% store indices of costs that were converted
if ~isempty(pwl1)
  om.userdata.pwl1 = pwl1;
end

% % [x, success, i] = nm.solve_power_flow(mpc, mpopt);
% [x, success, i] = nm.solve_opf(mpc, mpopt);
% opt = mpopt2nlpopt(mpopt, om.problem_type(), 'DEFAULT');
% [x, f, eflag, output, lambda] = om.solve(opt);
% keyboard;
