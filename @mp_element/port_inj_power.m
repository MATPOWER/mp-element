function [S, dS] = port_inj_power(mpe, x, idx)
%PORT_INJ_POWER  Evaluate port power injections and Jacobian.
%
%       S = mpe.port_inj_power(x, idx)
%       [S, dS] = mpe.port_inj_power(x, idx)
%       I = mpe.port_inj_current(x, idx)
%       [I, dI] = mpe.port_inj_current(x, idx)

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.
