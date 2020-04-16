function [F, J] = power_flow_equations(obj, x, va, vm, zr, zi, t)
%POWER_FLOW_EQUATIONS  Evaluates power flow equations
%   [F, J] = POWER_FLOW_EQUATIONS(OBJ, X, MPC)
%
%   Inputs
%       OBJ - 
%       X - 
%       VA - 
%       VM - 
%       Z - 
%       REF - 
%       PV - 
%       PQ - 
%
%   Returns
%       F - 
%       J -  
%
%   See also ...

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% index vector
pvq = [t.pv; t.pq];

%% update model state ([v_; z_]) from power flow state (x)
[v_, z_] = pfx2vz(obj, x, va, vm, zr, zi, t);

%% incidence matrix
C = obj.C;

%% Jacobian
if nargout > 1
    %% get port power injections with derivatives
    [S, Sva, Svm] = obj.port_inj_power([v_; z_], 1);

    SSva = C * Sva;
    SSvm = C * Svm;
    J = [   real(SSva(pvq,  pvq)) real(SSvm(pvq,  t.pq));
            imag(SSva(t.pq, pvq)) imag(SSvm(t.pq, t.pq))  ];
%     SSzr = C * Szr;
%     SSzi = C * Szi;
%     J = [
%         real(SSva(pvq,  pvq)) real(SSvm(pvq,  t.pq)) real(SSzr(pvq,  :)) real(SSzi(pvq,  :));
%         imag(SSva(t.pq, pvq)) imag(SSvm(t.pq, t.pq)) imag(SSzr(t.pq, :)) imag(SSzi(t.pq, :))  ];
else
    %% get port power injections (w/o derivatives)
    S = obj.port_inj_power([v_; z_], 1);
end

%% nodal power balance
SS = C * S;
F = [real(SS(pvq)); imag(SS(t.pq))];
