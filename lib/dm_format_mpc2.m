classdef dm_format_mpc2 < handle
%DM_FORMAT_MPC2  Mixin class for MATPOWER data model elements for MPC v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function tab = get_table(obj, dm)
            tab = dm.mpc.(obj.table);
        end
    end     %% methods
end         %% classdef
