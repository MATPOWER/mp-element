classdef mpx_reserves < mp_extension

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        function dmc_elements = dmc_element_classes(obj, dmc_class, fmt, mpopt)
            switch fmt
                case 'mpc2'
                    dmc_elements = { @dmce_reserve_gen_mpc2, ...
                                     @dmce_reserve_zone_mpc2 };
                otherwise
                    dmc_elements = {};      %% no modifications
            end
        end

        function dm_elements = dm_element_classes(obj, dm_class, task_tag, mpopt)
            switch task_tag
                case 'OPF'
                    dm_elements = { @dme_reserve_gen, ...
                                    @dme_reserve_zone };
                otherwise
                    dm_elements = {};       %% no modifications
            end
        end

        function nm_elements = nm_element_classes(obj, nm_class, task_tag, mpopt)
            switch task_tag
                case {'OPF'}
                    if strcmp(upper(mpopt.model), 'AC')
                        nm_elements = { @nme_reserve_gen_ac, ...
                                        @nme_reserve_zone_ac };
                    else
                        nm_elements = { @nme_reserve_gen_dc, ...
                                        @nme_reserve_zone_dc };
                    end
                otherwise
                    nm_elements = {};       %% no modifications
            end
        end

        function mm_elements = mm_element_classes(obj, mm_class, task_tag, mpopt)
            switch task_tag
                case {'OPF'}
                    mm_elements = { @mme_reserve_gen, ...
                                    @mme_reserve_zone };
                otherwise
                    mm_elements = {};       %% no modifications
            end
        end
    end     %% methods
end         %% classdef
