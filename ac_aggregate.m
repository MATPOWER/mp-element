classdef ac_aggregate < mp_aggregate% & ac_model
%AC_AGGREGATE Abstract class, explicitly a subclass of mp_aggregate and
%             implicitly assumed to be subclasses of ac_model as well

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
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

        function [G, Gv1, Gv2, Gzr, Gzi] = nodal_power_balance(obj, xx)
            %% node incidence matrix
            C = obj.getC();

            %% get port power injections with derivatives
            if nargout > 1
                [S, Sv1, Sv2, Szr, Szi] = obj.port_inj_power(xx, 1);
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

        function [g, dg] = opf_power_balance_fcn(obj, xc)
            xx = [xc{2} .* exp(1j*xc{1}); xc{3}+1j*xc{4}];
            if nargout > 1
                [G, Gv1, Gv2, Gzr, Gzi] = obj.nodal_power_balance(xx);
                dG = [Gv1 Gv2 Gzr Gzi];
                dg = [  real(dG);       %% P mismatch w.r.t v1, v2, zr, zi
                        imag(dG)    ];  %% Q mismatch w.r.t v1, v2, zr, zi
            else
                G = obj.nodal_power_balance(xx);
            end
            g = [ real(G);              %% active power (P) mismatch
                  imag(G) ];            %% reactive power (Q) mismatch
        end
    end     %% methods
end         %% classdef
