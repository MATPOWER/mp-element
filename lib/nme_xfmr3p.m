classdef nme_xfmr3p < nm_element & mp_form_acp

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
        function obj = nme_xfmr3p()
            obj@nm_element();
            obj.name = 'xfmr3p';
            obj.np = 6;             %% this is a 6 port element
        end

        function obj = build_params(obj, nm, dm)
            build_params@nm_element(obj, nm, dm);   %% call parent

            dme = obj.data_model_element(dm);
            bus_dme = dm.elements.(dme.cxn_type);
            nk = 3*obj.nk;

            base_kv = bus_dme.tab.base_kv(bus_dme.on(dme.fbus(dme.on)));
            z = dm.base_kva * (dme.r + 1j * dme.x) ./ dme.base_kva .* ...
                (dme.base_kv ./ ( base_kv / sqrt(3))) .^ 2;
            y = 1 ./ z;
            Y = [y;y;y];

            obj.Y = sparse( ...
                [1:nk 1:nk nk+1:2*nk nk+1:2*nk]', ...
                [1:nk nk+1:2*nk 1:nk nk+1:2*nk]', ...
                [Y; -Y; -Y; Y], 2*nk, 2*nk );
        end

        %%-----  PF methods  -----
        function obj = pf_data_model_update(obj, mm, nm, dm, mpopt)
            dme = obj.data_model_element(dm);
            pp = nm.get_idx('port');
            nn = obj.np/2;

            %% branch active power flow
            for p = 1:nn
                s_fr = nm.soln.gs_(pp.i1.xfmr3p(p):pp.iN.xfmr3p(p)) * dm.base_kva;
                s_to = nm.soln.gs_(pp.i1.xfmr3p(p+nn):pp.iN.xfmr3p(p+nn)) * dm.base_kva;

                %% update in the data model
                dme.tab.(sprintf('pl%d_fr', p))(dme.on) = real(s_fr);
                dme.tab.(sprintf('pf%d_fr', p))(dme.on) = cos(angle(s_fr));
                dme.tab.(sprintf('pl%d_to', p))(dme.on) = real(s_to);
                dme.tab.(sprintf('pf%d_to', p))(dme.on) = cos(angle(s_to));
            end
        end
    end     %% methods
end         %% classdef
