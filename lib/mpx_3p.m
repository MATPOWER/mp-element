classdef mpx_3p < mp_extension

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
                    dmc_elements = { ...
                        @mp.dmce_bus3p_mpc2, @mp.dmce_gen3p_mpc2, ...
                        @mp.dmce_load3p_mpc2, @mp.dmce_line3p_mpc2, ...
                        @mp.dmce_xfmr3p_mpc2, @mp.dmce_buslink_mpc2 ...
                    };
                otherwise
                    dmc_elements = {};
            end
        end

        function dm_elements = dm_element_classes(obj, dm_class, task_tag, mpopt)
            switch task_tag
                case {'PF', 'CPF'}
                    dm_elements = { ...
                        @dme_bus3p, @dme_gen3p, @dme_load3p, ...
                        @dme_line3p, @dme_xfmr3p, @dme_buslink ...
                    };
                case 'OPF'
                    dm_elements = { ...
                        @dme_bus3p_opf, @dme_gen3p_opf, @dme_load3p_opf, ...
                        @dme_line3p_opf, @dme_xfmr3p_opf, @dme_buslink_opf ...
                    };
                otherwise
                    dm_elements = {};
            end
        end

        function nm_elements = nm_element_classes(obj, nm_class, task_tag, mpopt)
            switch task_tag
                case {'PF', 'CPF'}
                    v_cartesian = mpopt.pf.v_cartesian;
                case {'OPF'}
                    v_cartesian = mpopt.opf.v_cartesian;
            end
            switch upper(mpopt.model)
                case 'AC'
                    if v_cartesian
                        nm_elements = { ...
                            @nme_bus3p_acc, @nme_gen3p_acc, @nme_load3p, ...
                            @nme_line3p, @nme_xfmr3p, @nme_buslink_acc ...
                        };
                    else
                        nm_elements = { ...
                            @nme_bus3p_acp, @nme_gen3p_acp, @nme_load3p, ...
                            @nme_line3p, @nme_xfmr3p, @nme_buslink_acp ...
                        };
                    end
                case 'DC'
                    nm_elements = {};       %% no modifications
            end
        end

        function mm_elements = mm_element_classes(obj, mm_class, task_tag, mpopt)
            switch task_tag
                case {'PF', 'CPF'}
                    switch upper(mpopt.model)
                        case 'AC'
                            if mpopt.pf.v_cartesian
                                mm_elements = { ...
                                    @mme_bus3p, @mme_gen3p, @mme_line3p, ...
                                    @mme_xfmr3p, @mme_buslink_pf_acc ...
                                };
                            else
                                mm_elements = { ...
                                    @mme_bus3p, @mme_gen3p, @mme_line3p, ...
                                    @mme_xfmr3p, @mme_buslink_pf_acp ...
                                };
                            end
                        case 'DC'
                            mm_elements = {};       %% no modifications
                    end
                case {'OPF'}
                    switch upper(mpopt.model)
                        case 'AC'
                            if mpopt.opf.v_cartesian
                                mm_elements = { ...
                                    @mme_bus3p_opf_acc, @mme_gen3p_opf, ...
                                    @mme_line3p_opf, @mme_xfmr3p_opf, ...
                                    @mme_buslink_opf_acc ...
                                };
                            else
                                mm_elements = { ...
                                    @mme_bus3p_opf_acp, @mme_gen3p_opf, ...
                                    @mme_line3p_opf, @mme_xfmr3p_opf, ...
                                    @mme_buslink_opf_acp ...
                                };
                            end
                        case 'DC'
                            mm_elements = {};       %% no modifications
                    end
            end
        end
    end     %% methods
end         %% classdef
