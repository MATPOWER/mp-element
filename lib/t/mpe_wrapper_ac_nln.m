classdef mpe_wrapper_ac_nln < handle

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpe = [];           %% wrapped mp_element object
    end
    
    methods
        function obj = mpe_wrapper_ac_nln_init(obj)
            obj.mpe = obj.mpe_class();      %% construct wrapped class
            
            %% copy base field values from wrapped object
            obj.name = obj.mpe.name;
            obj.mpc_field = obj.mpe.mpc_field;
            obj.np = obj.mpe.np;
            obj.nz = obj.mpe.nz;
        end
        
        function build_nln_params(obj, nm, mpc)
            %% build params for wrapped object
            obj.mpe.build_params(nm, mpc);
            
            %% remove other params
            obj.Y = [];
            obj.L = [];
            obj.M = [];
            obj.N = [];
            obj.i = [];
            obj.s = [];
            
            %% add nonlinear function/hessian params
            obj.inln = @(x_, sysx, idx)port_inj_current(obj.mpe, x_, sysx, idx);
            obj.snln = @(x_, sysx, idx)port_inj_power(obj.mpe, x_, sysx, idx);
            obj.inln_hess = @(x_, lam, sysx, idx)port_inj_current_hess(obj.mpe, x_, lam, sysx, idx);
            obj.snln_hess = @(x_, lam, sysx, idx)port_inj_power_hess(obj.mpe, x_, lam, sysx, idx);
        end

        function nk = count_nln(obj, mpc)
            obj.mpe.count(mpc);
        end
    end     %% methods
end         %% classdef
