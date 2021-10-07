classdef mme_pf_buslink_acc < mme_pf_buslink_ac

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'buslink';
%     end
    
    methods
        function obj = add_constraints(obj, mm, nm, dm, mpopt)
            nme = obj.network_model_element(nm);

            %% add constraints for matching
            %%  voltage angles at pv and pq nodes
            %%  voltage magnitudes at pq nodes
            [A_va_pq, A_va_pv, b_va, A_vm, b_vm] = pf_voltage_constraints(nme, mm.aux_data);

            %% prep variable set structs
            vs_va = struct('name', {'Vr_pq', 'Vr3_pq', 'Vr3_pq', 'Vr3_pq', ...
                                    'Vr_pv', 'Vr3_pv', 'Vr3_pv', 'Vr3_pv', ...
                                    'Vi_pq', 'Vi3_pq', 'Vi3_pq', 'Vi3_pq', ...
                                    'Vi_pv', 'Vi3_pv', 'Vi3_pv', 'Vi3_pv'}, ...
                            'idx', {{}, {1}, {2}, {3}, {}, {1}, {2}, {3}, ...
                                    {}, {1}, {2}, {3}, {}, {1}, {2}, {3}} );
            vs_vm = struct('name', {'Vr_pq', 'Vr3_pq', 'Vr3_pq', 'Vr3_pq', ...
                                    'Vi_pq', 'Vi3_pq', 'Vi3_pq', 'Vi3_pq'}, ...
                            'idx', {{}, {1}, {2}, {3}, {}, {1}, {2}, {3}} );

            fcn_va = @(xx)pf_va_fcn(nme, xx, [A_va_pq A_va_pv], b_va);
            mm.add_nln_constraint('buslink_va', length(b_va), 1, fcn_va, [], vs_va);

            fcn_vm = @(xx)pf_vm_fcn(nme, xx, A_vm, b_vm);
            mm.add_nln_constraint('buslink_vm', length(b_vm), 1, fcn_vm, [], vs_vm);

            %% call parent
            add_constraints@mme_pf_buslink_ac(obj, mm, nm, dm, mpopt);
        end
    end     %% methods
end         %% classdef
