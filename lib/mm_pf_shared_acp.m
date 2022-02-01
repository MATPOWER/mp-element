classdef mm_pf_shared_acp < mm_pf_shared_ac

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end
    
    methods
        function [vx_, z_, x_] = pf_convert_x(obj, mmx, nm, ad, only_v)
            %% x = obj.pf_convert(mmx, nm, ad)
            %% [v, z] = obj.pf_convert(mmx, nm, ad)
            %% [v, z, x] = obj.pf_convert(mmx, nm, ad)
            %% ... = obj.pf_convert(mmx, nm, ad, only_v)

            %% update v_, z_ from mmx
            nm_vars = nm.update_vars(mmx, ad);
            vx_ = nm_vars.vm .* exp(1j * nm_vars.va);
            z_ = nm_vars.zr + 1j * nm_vars.zi;

            %% update z, if requested
            if nargin < 5 || ~only_v
                z_ = obj.update_z(nm, vx_, z_, ad);
            end

            if nargout < 2
                vx_ = [vx_; z_];
            elseif nargout > 2
                x_ = [vx_; z_];
            end
        end
    end     %% methods
end         %% classdef
