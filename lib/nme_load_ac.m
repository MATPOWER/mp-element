classdef nme_load_ac < nme_load% & mp_form_ac

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function obj = build_params(obj, nm, dm)
            build_params@nme_load(obj, nm, dm);     %% call parent

            dme = obj.data_model_element(dm);
            obj.s = dme.Pd + 1j * dme.Qd;           %% complex power demand

            %% experimental system-wide ZIP loads (for backward compatibility)
            [s, Sd, Y] = dme.sys_wide_zip_loads(dm);
            if ~isempty(s)
                obj.s = s;
                obj.Y = Y;
%                 obj.inln = @(x_, sysx, idx)port_inj_current_nln(obj, Sd, x_, sysx, idx);
                obj.snln = @(x_, sysx, idx)port_inj_power_nln(obj, Sd, x_, sysx, idx);
            end
        end

%         function [I, Iv1, Iv2, Izr, Izi] = port_inj_current_nln(obj, Sd, x_, sysx, idx)
%             if nargin < 5
%                 idx = [];
%                 if nargin < 4
%                     sysx = 1;
%                 end
%             end
% 
%             [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);
%             if isempty(idx)
%                 Sdi = Sd;
%             else
%                 Sdi = Sd(idx);
%             end
%             S = abs(vi_) .* Sdi;
%             I = S ./ vi_;
% 
%             if nargout > 1
%                 nv = length(v_);
%                 nz = length(z_);
%                 ni = length(S);
%                 if isempty(idx)
%                     idx = (1:ni);
%                 end
%                 Sv1 = sparse(ni, nv);
%                 Sv2 = sparse(1:ni, idx, Sdi, ni, nv);
%                 if nargout > 3
%                     Szr = sparse(ni, nz);
%                     Szi = Szr;
%                 end
%                 if sysx
%                     Ct = obj.C';
%                     Sv1 = Sv1 * Ct;
%                     Sv2 = Sv2 * Ct;
%                     if nargout > 3  %% Szr, Szi are empty, but num of rows is needed
%                         Dt = obj.D';
%                         Szr = Szr * Dt;
%                         Szi = Szi * Dt;
%                     end
%                 end
%             end
%         end

        function [S, Sv1, Sv2, Szr, Szi] = port_inj_power_nln(obj, Sd, x_, sysx, idx)
            if nargin < 5
                idx = [];
                if nargin < 4
                    sysx = 1;
                end
            end

            [v_, z_, vi_] = obj.x2vz(x_, sysx, idx);
            if isempty(idx)
                Sdi = Sd;
            else
                Sdi = Sd(idx);
            end
            S = abs(vi_) .* Sdi;

            if nargout > 1
                nv = length(v_);
                nz = length(z_);
                ni = length(S);
                if isempty(idx)
                    idx = (1:ni);
                end
                Sv1 = sparse(ni, nv);
                Sv2 = sparse(1:ni, idx, Sdi, ni, nv);
                if nargout > 3
                    Szr = sparse(ni, nz);
                    Szi = Szr;
                end
                if sysx
                    Ct = obj.C';
                    Sv1 = Sv1 * Ct;
                    Sv2 = Sv2 * Ct;
                    if nargout > 3  %% Szr, Szi are empty, but num of rows is needed
                        Dt = obj.D';
                        Szr = Szr * Dt;
                        Szi = Szi * Dt;
                    end
                end
            end
        end

        %%-----  CPF methods  -----
        function obj = cpf_data_model_update(obj, mm, nm, dm, mpopt)
            ad = mm.get_userdata('aux_data');
            dme = obj.data_model_element(dm);
            dm = dme.parameterized(dm, ad.dmb, ad.dmt, mm.soln.x(end));
        end
    end     %% methods
end         %% classdef
