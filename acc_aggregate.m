classdef acc_aggregate < ac_aggregate & mp_model_acc

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        vr = [];
        vi = [];
    end
    
    methods
        function obj = acc_aggregate()
            obj@ac_aggregate();
            obj.element_classes = ...
                { @acc_bus, @acc_gen, @acc_load, @acc_branch, @acc_shunt };
            if isempty(obj.node)    %% skip if constructed from existing object
                obj.init_set_types();   %% should be called in mp_idx_manager
                                        %% constructor, if not for:
            end                         %% https://savannah.gnu.org/bugs/?52614
        end

        function obj = def_set_types(obj)
            def_set_types@ac_aggregate(obj);        %% call parent first
            obj.set_types.vr = 'REAL VOLTAGE VARS (vr)';
            obj.set_types.vi = 'IMAG VOLTAGE VARS (vi)';
        end


        %%-----  OPF methods  -----
        function x_ = x2x_(obj, x)
            %% convert (real) opt_model x to (complex) network model x_
            nv_ = obj.nv / 2;       %% number of voltage vars (sysx=1)
            nz_ = obj.nz;           %% number of state vars
            vr = x(1:nv_, :);       b = nv_;
            vi = x(b+1:b+nv_, :);   b = b + nv_;
            zr = x(b+1:b+nz_, :);   b = b + nz_;
            zi = x(b+1:b+nz_, :);
            x_ = [vr+1j*vi; zr+1j*zi];
        end

        function names = opf_legacy_user_var_names(obj)
            names = {'Vr', 'Vi', 'Pg', 'Qg'};
        end
    end     %% methods
end         %% classdef
