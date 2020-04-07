classdef ac_aggregate < mp_aggregate% & ac_model
%AC_AGGREGATE Abstract class, explicitly a subclass of mp_aggregate and
%             implicitly assumed to be subclasses of ac_model as well

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        zr = [];
        zi = [];
    end
    
    methods
        %% constructor
        function obj = ac_aggregate(varargin)
            obj@mp_aggregate(varargin{:});
        end

        function obj = def_set_types(obj)
            def_set_types@mp_aggregate(obj);        %% call parent first
            obj.set_types.zr = 'NON-VOLTAGE VARS REAL (zr)';
            obj.set_types.zi = 'NON-VOLTAGE VARS IMAG (zi)';
        end

        function obj = build_params(obj, asm, mpc)
            %% call parent to build individual element parameters
            build_params@mp_aggregate(obj, asm, mpc);

            %% aggregate parameters from individual elements
            obj.Y = obj.stack_matrix_params('Y', 1);
            obj.L = obj.stack_matrix_params('L', 0);
            obj.M = obj.stack_matrix_params('M', 1);
            obj.N = obj.stack_matrix_params('N', 0);
            obj.i = obj.stack_vector_params('i');
            obj.s = obj.stack_vector_params('s');
        end

        function [G, Gv1, Gv2, Gzr, Gzi] = nodal_complex_current_balance(obj, x_)
            %% node incidence matrix
            C = obj.getC();

            %% get port current injections with derivatives
            if nargout > 1
                [I, Iv1, Iv2, Izr, Izi] = obj.port_inj_current(x_, 1);
                Gv1 = C * Iv1;      %% Gva or Gvr
                Gv2 = C * Iv2;      %% Gvm or Gvi
                Gzr = C * Izr;
                Gzi = C * Izi;
            else
                I = obj.port_inj_current(x, 1);
            end

            %% nodal current balance
            G = C * I;
        end

        function [G, Gv1, Gv2, Gzr, Gzi] = nodal_complex_power_balance(obj, x_)
            %% node incidence matrix
            C = obj.getC();

            %% get port power injections with derivatives
            if nargout > 1
                [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(x_, 1);
                Gv1 = C * Sv1;      %% Gva or Gvr
                Gv2 = C * Sv2;      %% Gvm or Gvi
                Gzr = C * Szr;
                Gzi = C * Szi;
            else
                S = obj.port_inj_power(x, 1);
            end

            %% nodal power balance
            G = C * S;
        end

        function d2G = nodal_complex_current_balance_hess(obj, x_, lam)
            %% node incidence matrix
            Ct = obj.getC('tr');

            %% get port power injection hessians
            d2G = obj.port_inj_current_hess(x_, Ct * lam);
        end

        function d2G = nodal_complex_power_balance_hess(obj, x_, lam)
            %% node incidence matrix
            Ct = obj.getC('tr');

            %% get port power injection hessians
            d2G = obj.port_inj_power_hess(x_, Ct * lam);
        end

        function [g, dg] = opf_current_balance_fcn(obj, x)
            x_ = obj.x2x_(x);           %% convert real to complex x
            if nargout > 1
                [G, Gv1, Gv2, Gzr, Gzi] = obj.nodal_complex_current_balance(x_);
                Gx = [Gv1 Gv2 Gzr Gzi];
                dg = [  real(Gx);       %% Re{I} mismatch w.r.t v1, v2, zr, zi
                        imag(Gx)    ];  %% Im{I} mismatch w.r.t v1, v2, zr, zi
            else
                G = obj.nodal_complex_current_balance(x_);
            end
            g = [ real(G);              %% real current mismatch
                  imag(G) ];            %% imaginary current mismatch
        end

        function [g, dg] = opf_power_balance_fcn(obj, x)
            x_ = obj.x2x_(x);           %% convert real to complex x
            if nargout > 1
                [G, Gv1, Gv2, Gzr, Gzi] = obj.nodal_complex_power_balance(x_);
                Gx = [Gv1 Gv2 Gzr Gzi];
                dg = [  real(Gx);       %% P mismatch w.r.t v1, v2, zr, zi
                        imag(Gx)    ];  %% Q mismatch w.r.t v1, v2, zr, zi
            else
                G = obj.nodal_complex_power_balance(x_);
            end
            g = [ real(G);              %% active power (P) mismatch
                  imag(G) ];            %% reactive power (Q) mismatch
        end

        function d2G = opf_current_balance_hess(obj, x, lam)
            x_ = obj.x2x_(x);           %% convert real to complex x
            nlam = length(lam) / 2;
            lamIr = lam(1:nlam);
            lamIi = lam((1:nlam)+nlam);

            d2Gr = obj.nodal_complex_current_balance_hess(x_, lamIr);
            d2Gi = obj.nodal_complex_current_balance_hess(x_, lamIi);

            d2G = real(d2Gr) + imag(d2Gi);
        end

        function d2G = opf_power_balance_hess(obj, x, lam)
            x_ = obj.x2x_(x);           %% convert real to complex x
            nlam = length(lam) / 2;
            lamP = lam(1:nlam);
            lamQ = lam((1:nlam)+nlam);

            d2Gr = obj.nodal_complex_power_balance_hess(x_, lamP);
            d2Gi = obj.nodal_complex_power_balance_hess(x_, lamQ);

            d2G = real(d2Gr) + imag(d2Gi);
        end
    end     %% methods
end         %% classdef
