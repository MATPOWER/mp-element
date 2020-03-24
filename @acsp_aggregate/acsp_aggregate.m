classdef acsp_aggregate < ac_aggregate & acsp_model

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
        %% constructor
        function obj = acsp_aggregate(varargin)
            obj@ac_aggregate(varargin{:});
            obj.element_classes = ...
                { @acsp_bus, @ac_gen, @ac_load, @ac_branch };
            if isempty(obj.node)    %% skip if constructed from existing object
                obj.init_set_types();   %% should be called in mp_idx_manager
                                        %% constructor, if not for:
                                        %% https://savannah.gnu.org/bugs/?52614
            end
        end

        function obj = def_set_types(obj)
            def_set_types@ac_aggregate(obj);        %% call parent first
            obj.set_types.va = 'VOLTAGE ANG VARS (va)';
            obj.set_types.vm = 'VOLTAGE MAG VARS (vm)';
        end

        function x_ = x2x_(obj, x)
            %% convert real to complex x
            nx_ = size(x, 1) / 2;   %% number of state vars
            nv_ = obj.nv / 2;       %% number of voltage vars (sysx=1)
            nz_ = nx_ - nv_;        %% number of non-voltage state vars
            va = x(1:nv_, :);       b = nv_;
            vm = x(b+1:b+nv_, :);   b = b + nv_;
            zr = x(b+1:b+nz_, :);   b = b + nz_;
            zi = x(b+1:b+nz_, :);
            x_ = [vm .* exp(1j*va); zr+1j*zi];
        end

        function d2G = nodal_power_balance_hess(obj, x_, lam)
            %% node incidence matrix
            Ct = obj.getC('tr');

            %% get port power injection hessians
            d2G = obj.port_inj_power_hess(x_, Ct * lam);
        end

        function d2G = opf_power_balance_hess(obj, x, lam)
            x_ = obj.x2x_(x);           %% convert real to complex x
            nlam = length(lam) / 2;
            lamP = lam(1:nlam);
            lamQ = lam((1:nlam)+nlam);

            d2Gr = obj.nodal_power_balance_hess(x_, lamP);
            d2Gi = obj.nodal_power_balance_hess(x_, lamQ);

            d2G = real(d2Gr) + imag(d2Gi);
        end

        function add_opf_node_balance_constraints(obj, om)
            %% power balance constraints
            nn = obj.node.N;            %% number of nodes
            fcn_mis = @(x)opf_power_balance_fcn(obj, x);
            hess_mis = @(x, lam)opf_power_balance_hess(obj, x, lam);
            om.add_nln_constraint({'Pmis', 'Qmis'}, [nn;nn], 1, fcn_mis, hess_mis);
        end
    end     %% methods
end         %% classdef
