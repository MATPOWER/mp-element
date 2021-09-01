classdef mp_network_acps_node_test < mp_network_acps

%   MATPOWER
%   Copyright (c) 2019-2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    methods
        %% constructor
        function obj = mp_network_acps_node_test()
            obj@mp_network_acps();
            obj.element_classes = ...
                { @nme_bus_nld_acp_node_test, @nme_bus_ld_acp_node_test, ...
                    @nme_gen_acp, @nme_branch_acp };

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  INIT_SET_TYPES() is called explicitly as needed
            %%              (if obj.node is empty) in BUILD() and DISPLAY(),
            %%              after object construction, but before object use.
        end

        %%-----  OPF methods  -----
        function opt = opf_solve_opts(obj, mm, dm, mpopt)
            opt = mpopt2nlpopt(mpopt, mm.problem_type());

            if mpopt.opf.start < 2
                %% initialize interior point
                x0 = obj.opf_interior_x0(mm, dm);

                %% set voltages
                %% va equal to angle of 1st ref bus
                %% vm equal to avg of clipped limits
                vv = mm.get_idx();
                gen_dme = dm.elements.gen;
                varefs = [];
                for k = gen_dme.nbet:-1:1
                    if dm.elements.is_index_name(gen_dme.cxn_type{k})
                        bus_dme{k} = dm.elements.(gen_dme.cxn_type{k});
                        varefs_k = bus_dme{k}.va_start(find(bus_dme{k}.type == NODE_TYPE.REF));
                        varefs = [varefs_k; varefs];
                    else
                        bus_dme{k} = [];
                    end
                end
                for k = 1:gen_dme.nbet
                    if ~isempty(bus_dme{k})
                        vm_ub = min(bus_dme{k}.vm_ub, 1.5);
                        vm_lb = max(bus_dme{k}.vm_lb, 0.5);
                        vm = (vm_ub + vm_lb) / 2;
                        vVa = ['va_' gen_dme.cxn_type{k}];
                        vVm = ['vm_' gen_dme.cxn_type{k}];
                        x0(vv.i1.(vVa):vv.iN.(vVa)) = varefs(1);%% angles set to first reference angle
                        x0(vv.i1.(vVm):vv.iN.(vVm)) = vm;       %% voltage magnitudes
                    end
                end

                opt.x0 = x0;
            end
        end
    end     %% methods
end         %% classdef
