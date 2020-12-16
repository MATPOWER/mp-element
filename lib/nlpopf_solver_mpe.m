function [results, success, raw] = nlpopf_solver_mpe(opf, mpopt)

%%----- initialization -----
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% unpack data
om = opf.mm;
mpc = om.get_mpc();
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
[vv, ll, nne, nni] = om.get_idx();

%% problem dimensions
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ny = om.getN('var', 'y');   %% number of piece-wise linear costs

[x, f, eflag, output, lambda, success] = deal(om.soln.x, om.soln.f, ...
    om.soln.eflag, om.soln.output, om.soln.lambda, opf.success);

%% find branches with flow limits
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
nl2 = length(il);           %% number of constrained lines

%% update solution data
if mpopt.opf.v_cartesian
    Vi = x(vv.i1.Vi:vv.iN.Vi);
    Vr = x(vv.i1.Vr:vv.iN.Vr);
    V = Vr + 1j*Vi;
    Va = angle(V);
    Vm = abs(V);
else
    Va = x(vv.i1.Va:vv.iN.Va);
    Vm = x(vv.i1.Vm:vv.iN.Vm);
    V = Vm .* exp(1j*Va);
end
Pg = x(vv.i1.Pg:vv.iN.Pg);
Qg = x(vv.i1.Qg:vv.iN.Qg);

%%-----  calculate return values  -----
%% update voltages & generator outputs
bus(:, VA) = Va * 180/pi;
bus(:, VM) = Vm;
gen(:, PG) = Pg * baseMVA;
gen(:, QG) = Qg * baseMVA;
gen(:, VG) = Vm(gen(:, GEN_BUS));

%% compute branch flows
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);    %% build admittance matrices
Sf = V(branch(:, F_BUS)) .* conj(Yf * V);   %% cplx pwr at "from" bus, p.u.
St = V(branch(:, T_BUS)) .* conj(Yt * V);   %% cplx pwr at "to" bus, p.u.
branch(:, PF) = real(Sf) * baseMVA;
branch(:, QF) = imag(Sf) * baseMVA;
branch(:, PT) = real(St) * baseMVA;
branch(:, QT) = imag(St) * baseMVA;

%% line constraint is typically on square of limit
%% so we must fix multipliers
muSf = zeros(nl, 1);
muSt = zeros(nl, 1);
if ~isempty(il)
    if upper(mpopt.opf.flow_lim(1)) == 'P'
        muSf(il) = lambda.ineqnonlin(nni.i1.Sf:nni.iN.Sf);
        muSt(il) = lambda.ineqnonlin(nni.i1.St:nni.iN.St);
    else
        muSf(il) = 2 * lambda.ineqnonlin(nni.i1.Sf:nni.iN.Sf) .* branch(il, RATE_A) / baseMVA;
        muSt(il) = 2 * lambda.ineqnonlin(nni.i1.St:nni.iN.St) .* branch(il, RATE_A) / baseMVA;
    end
end

%% update Lagrange multipliers
if mpopt.opf.v_cartesian
    if om.userdata.veq
        lam = lambda.eqnonlin(nne.i1.Veq:nne.iN.Veq);
        mu_Vmax = zeros(size(lam));
        mu_Vmin = zeros(size(lam));
        mu_Vmax(lam > 0) =  lam(lam > 0);
        mu_Vmin(lam < 0) = -lam(lam < 0);
        bus(om.userdata.veq, MU_VMAX) = mu_Vmax;
        bus(om.userdata.veq, MU_VMIN) = mu_Vmin;
    end
    bus(om.userdata.viq, MU_VMAX) = lambda.ineqnonlin(nni.i1.Vmax:nni.iN.Vmax);
    bus(om.userdata.viq, MU_VMIN) = lambda.ineqnonlin(nni.i1.Vmin:nni.iN.Vmin);
else
    bus(:, MU_VMAX)  = lambda.upper(vv.i1.Vm:vv.iN.Vm);
    bus(:, MU_VMIN)  = lambda.lower(vv.i1.Vm:vv.iN.Vm);
end
gen(:, MU_PMAX)  = lambda.upper(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_PMIN)  = lambda.lower(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_QMAX)  = lambda.upper(vv.i1.Qg:vv.iN.Qg) / baseMVA;
gen(:, MU_QMIN)  = lambda.lower(vv.i1.Qg:vv.iN.Qg) / baseMVA;
if mpopt.opf.current_balance
    %% convert current balance shadow prices to equivalent lamP and lamQ
    %% P + jQ = (Vr + jVi) * (M - jN)
    %% M = (Vr P + Vi Q) / (Vr^2 + Vi^2)
    %% N = (Vi P - Vr Q) / (Vr^2 + Vi^2)
    %% lamP = df/dP = df/dM * dM/dP + df/dN + dN/dP
    %% lamQ = df/dQ = df/dM * dM/dQ + df/dN + dN/dQ
    VV = V ./ (V .* conj(V));   %% V / Vm^2
    VVr = real(VV);
    VVi = imag(VV);
    lamM = lambda.eqnonlin(nne.i1.rImis:nne.iN.rImis);
    lamN = lambda.eqnonlin(nne.i1.iImis:nne.iN.iImis);
    bus(:, LAM_P) = (VVr.*lamM + VVi.*lamN) / baseMVA;
    bus(:, LAM_Q) = (VVi.*lamM - VVr.*lamN) / baseMVA;
else
    bus(:, LAM_P) = lambda.eqnonlin(nne.i1.Pmis:nne.iN.Pmis) / baseMVA;
    bus(:, LAM_Q) = lambda.eqnonlin(nne.i1.Qmis:nne.iN.Qmis) / baseMVA;
end
branch(:, MU_SF) = muSf / baseMVA;
branch(:, MU_ST) = muSt / baseMVA;

%% package up results
nlnN = 2*nb + 2*nl;     %% because muSf and muSt are nl x 1, not nl2 x 1

%% extract multipliers for nonlinear constraints
kl = find(lambda.eqnonlin(1:2*nb) < 0);
ku = find(lambda.eqnonlin(1:2*nb) > 0);
nl_mu_l = zeros(nlnN, 1);
nl_mu_u = [zeros(2*nb, 1); muSf; muSt];
nl_mu_l(kl) = -lambda.eqnonlin(kl);
nl_mu_u(ku) =  lambda.eqnonlin(ku);

if isfield(lambda, 'ineqnonlin')
    lam_nli = lambda.ineqnonlin;
else
    lam_nli = [];
end

mu = struct( ...
  'var', struct('l', lambda.lower, 'u', lambda.upper), ...
  'nln', struct('l', nl_mu_l, 'u', nl_mu_u), ...
  'nle', lambda.eqnonlin, ...
  'nli', lam_nli, ...
  'lin', struct('l', lambda.mu_l, 'u', lambda.mu_u) );

results = mpc;
[results.bus, results.branch, results.gen, ...
    results.om, results.x, results.mu, results.f] = ...
        deal(bus, branch, gen, om, x, mu, f);

pimul = [ ...
  results.mu.nln.l - results.mu.nln.u;
  results.mu.lin.l - results.mu.lin.u;
  -ones(ny>0, 1);
  results.mu.var.l - results.mu.var.u;
];
raw = struct('xr', x, 'pimul', pimul, 'info', eflag, 'output', output);
