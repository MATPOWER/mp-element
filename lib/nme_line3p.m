classdef nme_line3p < nm_element & mp_form_acp

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'branch';
%     end
    
    methods
        %% constructor
        function obj = nme_line3p()
            obj@nm_element();
            obj.name = 'line3p';
            obj.np = 6;             %% this is a 6 port element
        end

        function obj = build_params(obj, nm, dm)
            build_params@nm_element(obj, nm, dm);   %% call parent

            dme = obj.data_model_element(dm);
            bus_dme = dm.elements.(dme.cxn_type);
            nk = 3*obj.nk;

            base_kv = bus_dme.tab.base_kv(bus_dme.on(dme.fbus(dme.on))) / sqrt(3);
            base_z = 1000 / dm.base_kva * base_kv .^ 2;

            sf = (base_z ./ dme.len) * ones(1, 6);      %% scale factor
            y_od = dme.ys(dme.lc, :) .* sf;             %% - of off-diagonal
            y_d  = y_od + dme.yc(dme.lc, :)/2 .* sf;    %% diagonal

            Y_od = obj.vec2symmat_stacked(y_od);        %% - of off-diagonal
            Y_d  = obj.vec2symmat_stacked(y_d);         %% diagonal

            obj.Y = [ Y_d   -Y_od;
                     -Y_od   Y_d ];
        end

        function M = vec2symmat_stacked(obj, vv)
            %% analogous to making a symmetric matrix from a vector of 6 values
            %% M = [v(1) v(2) v(3);
            %%      v(2) v(4) v(5);
            %%      v(3) v(5) v(6) ];
            %% except v(k) is now a diagonal matrix build from vector vv(:, k)
            n = size(vv, 1);
            v = cell(6, 1);
            for k = 1:6
                v{k} = spdiags(vv(:, k), 0, n, n);
            end
            M = [v{1} v{2} v{3};
                 v{2} v{4} v{5};
                 v{3} v{5} v{6} ];
        end

        %%-----  PF methods  -----
        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            dme = obj.data_model_element(dm);
            pp = nm.get_idx('port');
            nn = obj.np/2;

            %% branch active power flow
            for p = 1:nn
                s_fr = nm.soln.gs_(pp.i1.line3p(p):pp.iN.line3p(p)) * dm.base_kva;
                s_to = nm.soln.gs_(pp.i1.line3p(p+nn):pp.iN.line3p(p+nn)) * dm.base_kva;

                %% update in the data model
                dme.tab.(sprintf('pl%d_fr', p))(dme.on) = real(s_fr);
                dme.tab.(sprintf('pf%d_fr', p))(dme.on) = cos(angle(s_fr));
                dme.tab.(sprintf('pl%d_to', p))(dme.on) = real(s_to);
                dme.tab.(sprintf('pf%d_to', p))(dme.on) = cos(angle(s_to));
            end
        end
    end     %% methods
end         %% classdef
