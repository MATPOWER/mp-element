classdef nme_buslink_acp < nme_buslink & mp_form_acp

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end

    methods
        %%-----  PF methods  -----
        function obj = pf_add_constraints(obj, mm, nm, dm, mpopt)
            %% add constraints for matching
            %%  voltage angles at pv and pq nodes
            %%  voltage magnitudes at pq nodes
            [A_va, b_va, A_vm, b_vm] = pf_voltage_constraints(obj, mm.aux_data);

            %% prep variable set structs
            vs_va = struct('name', {'Va_pv', 'Va3_pv', 'Va3_pv', 'Va3_pv', ...
                                    'Va_pq', 'Va3_pq', 'Va3_pq', 'Va3_pq'}, ...
                            'idx', {{}, {1}, {2}, {3}, {}, {1}, {2}, {3}} );
            vs_vm = struct('name', {'Vm_pq', 'Vm3_pq', 'Vm3_pq', 'Vm3_pq'}, ...
                            'idx', {{}, {1}, {2}, {3}} );

            %% add constraints
            mm.add_lin_constraint('buslink_va', A_va, b_va, b_va, vs_va);
            mm.add_lin_constraint('buslink_vm', A_vm, b_vm, b_vm, vs_vm);
        end
    end     %% methods
end         %% classdef
