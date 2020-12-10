function om = setup_opf(obj, dm, mpopt)
%SETUP_OPF  Build the OPF optimization model for this network model
%   OM = OBJ.SETUP_OPF(MPC, MPOPT)
%
%   Inputs
%       OBJ - 
%       MPC - 
%       MPOPT - 
%
%   Returns
%       OM - 
%
%   See also ...

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% MATPOWER options
if nargin < 3
    mpopt = mpoption;
end

%% create optimization model
% om = opt_model();
mpc = dm.mpc;
om = opf_model(mpc);    %% switch back to simple opt_model, if possible
if obj.np ~= 0      %% skip for empty model
    obj.add_opf_vars(obj, om, dm, mpopt);
    obj.add_opf_constraints(obj, om, dm, mpopt);
    obj.add_opf_costs(obj, om, dm, mpopt);
    
    %% execute userfcn callbacks for 'formulation' stage
    if isfield(mpc, 'userfcn')
        userfcn = mpc.userfcn;
    else
        userfcn = [];
    end
    om = run_userfcn(userfcn, 'formulation', om, mpopt);
end
