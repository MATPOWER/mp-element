classdef mp_network_dc < mp_network & mp_form_dc

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        va = [];
        z  = [];
    end
    
    methods
        %% constructor
        function obj = mp_network_dc()
            obj@mp_network();
            obj.element_classes = { ...
                @nme_bus_dc, @nme_gen_dc, @nme_load_dc, ...
                    @nme_branch_dc, @nme_shunt_dc };

            %% Due to a bug related to inheritance in constructors in
            %% Octave 5.2 and earlier (https://savannah.gnu.org/bugs/?52614),
            %% INIT_SET_TYPES() cannot be called directly in the
            %% MP_IDX_MANAGER constructor, as desired.
            %%
            %% WORKAROUND:  INIT_SET_TYPES() is called explicitly as needed
            %%              (if obj.node is empty) in BUILD() and DISPLAY(),
            %%              after object construction, but before object use.
        end

        function obj = def_set_types(obj)
            def_set_types@mp_network(obj);      %% call parent first
            obj.set_types.va = 'VOLTAGE VARS (va)';
            obj.set_types.z  = 'NON-VOLTAGE VARS (z)';
        end

        function obj = build_params(obj, nm, dm)
            %% call parent to build individual element parameters
            build_params@mp_network(obj, nm, dm);

            %% aggregate parameters from individual elements
            obj.B = obj.stack_matrix_params('B', 1);
            obj.K = obj.stack_matrix_params('K', 0);
            obj.p = obj.stack_vector_params('p');
        end

        function obj = port_inj_soln(obj)
            %% compute port injections
            obj.soln.gp = obj.port_inj_power(obj.soln.x);
        end


        %%-----  PF methods  -----
        function [vx, z, x] = pf_convert_x(obj, mmx, ad, only_v)
            %% x = obj.pf_convert(mmx, ad)
            %% [v, z] = obj.pf_convert(mmx, ad)
            %% [v, z, x] = obj.pf_convert(mmx, ad)
            %% ... = obj.pf_convert(mmx, ad, only_v)

            %% update v_, z_ from mmx
            vx = ad.va;
            vx([ad.pv; ad.pq]) = mmx(1:ad.npv+ad.npq);      %% va
            z = ad.z;

            %% update z, if requested
            if nargin < 4 || ~only_v
                z = obj.pf_update_z(vx, z, ad);
            end

            %% prepare return values
            if nargout < 2
                vx = [vx; z];
            elseif nargout > 2
                x = [vx; z];
            end
        end

        function z = pf_update_z(obj, v, z, ad)
            %% update/allocate slack node active power injections
            
            %% coefficient matrix for power injection states at slack bus
            CC = obj.C(ad.ref, :) * obj.get_params([], 'K') * obj.D';
            jr = find(any(CC, 1));  %% indices of corresponding states

            %% power injections at slack nodes
            idx = find(any(obj.C(ad.ref, :), 1));  %% ports connected to slack nodes
            Pref = obj.C(ad.ref, idx) * obj.port_inj_power([v; z], 1, idx);

            %% allocate active power at slack nodes to 1st direct inj state
            %% find all z (except first one) with direct injection at each
            %% slack node
            [i, j] = find(CC);
            if size(i, 2) > 1, i = i'; j = j'; end
            ij = sortrows([i j]);       %% 1st state comes 1st for each node
            [~, k1] = unique(ij(:, 1), 'first');%% index of 1st entry for each node
            %% all included states that are not 1st at their node
            jn = unique(ij(~ismember(1:length(i), k1), 2));

            %% if we have extra states (more than 1) for any node(s)
            if ~isempty(jn)
                %% augment update equation CC * (z - zprev) = -Pref with
                %% additional rows to force these states to remain fixed
                I = speye(obj.nz);
                CC = [CC; I(jn, :)];
                Pref = [Pref; zeros(length(jn), 1)];
            end

            %% update z for active injections at slack nodes
            z(jr) = z(jr) - CC(:, jr) \ Pref;
        end


        %%-----  OPF methods  -----
        function [vx, z, x] = opf_convert_x(obj, mmx, ad)
            nm_vars = obj.update_vars(mmx, ad);

            %% convert (real) math model x to network model x
            vx = nm_vars.va;
            z  = nm_vars.z;
            if nargout < 2
                vx = [vx; z];
            elseif nargout > 2
                x = [vx; z];
            end
        end

        function names = opf_legacy_user_var_names(obj)
            names = {'Va', 'Pg'};
        end
    end     %% methods
end         %% classdef
