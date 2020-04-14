classdef acp_aggregate < ac_aggregate & acp_model

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        va = [];
        vm = [];
    end
    
    methods
        function obj = def_set_types(obj)
            def_set_types@ac_aggregate(obj);        %% call parent first
            obj.set_types.va = 'VOLTAGE ANG VARS (va)';
            obj.set_types.vm = 'VOLTAGE MAG VARS (vm)';
        end

        function x_ = x2x_(obj, x)
            %% convert (real) opt_model x to (complex) network model x_
            nv_ = obj.nv / 2;       %% number of voltage vars (sysx=1)
            nz_ = obj.nz;           %% number of state vars
            va = x(1:nv_, :);       b = nv_;
            vm = x(b+1:b+nv_, :);   b = b + nv_;
            zr = x(b+1:b+nz_, :);   b = b + nz_;
            zi = x(b+1:b+nz_, :);
            x_ = [vm .* exp(1j*va); zr+1j*zi];
        end
    end     %% methods
end         %% classdef