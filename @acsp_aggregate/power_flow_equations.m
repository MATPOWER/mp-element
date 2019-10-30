function [F, J] = power_flow_equations(obj, x, va, vm, z, ref, pv, pq)
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

npv = length(pv);
npq = length(pq);
nv = length(va);
nz = length(z);
pvq = [pv; pq];

%% update model state ([v; z]) from power flow state (x)
va(pvq) = x(1:npv+npq);
vm(pq) = x(npv+npq+1:end);
v = vm .* exp(1j * va);

%% get port power injections with derivatives
[S, Sva, Svm, Szr, Szi] = port_inj_power(obj, [v; z], 1);

%% nodal power balance
C = obj.getC();
SS = C * S;
F = [real(SS(pvq)); imag(SS(pq))];

%% Jacobian
if nargout > 1
    SSva = C * Sva;
    SSvm = C * Svm;
    J = [   real(SSva(pvq, pvq)) real(SSvm(pvq, pq));
            imag(SSva(pq,  pvq)) imag(SSvm(pq,  pq))  ];
%     SSzr = C * Szr;
%     SSzi = C * Szi;
%     J = [
%         real(SSva(pvq, pvq)) real(SSvm(pvq, pq)) real(SSzr(pvq, :)) real(SSzi(pvq, :));
%         imag(SSva(pq,  pvq)) imag(SSvm(pq,  pq)) imag(SSzr(pq,  :)) imag(SSzi(pq,  :))  ];
end