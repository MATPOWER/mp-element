classdef mp_network_acp < mp_network_ac% & mp_form_acp

%   MATPOWER
%   Copyright (c) 2019-2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        va = [];
        vm = [];
    end
    
    methods
        function obj = mp_network_acp()
            obj@mp_network_ac();
            obj.element_classes = ...
                { @nme_bus_acp, @nme_gen_acp, @nme_load_acp, ...
                    @nme_branch_acp, @nme_shunt_acp };

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
            def_set_types@mp_network_ac(obj);   %% call parent first
            obj.set_types.va = 'VOLTAGE ANG VARS (va)';
            obj.set_types.vm = 'VOLTAGE MAG VARS (vm)';
        end

        function va = initial_voltage_angle(obj, idx)
            va = obj.params_var('va');
            if nargin > 1 && ~isempty(idx)
                va = va(idx);
            end
        end

        function [va, vm] = va_vm(obj, v1, v2)
            va = v1;
            vm = v2;
        end

        %%-----  CPF methods  -----
        function efv = cpf_event_vlim(obj, cx, opt, mm, dm, mpopt)
            %% convert cx.x back to v_
            ad = mm.get_userdata('aux_data');

            %% get current node voltage magnitudes and bounds
            [v_, ~] = obj.cpf_convert_x(cx.x, ad, 1);
            [~, vm_min, vm_max] = obj.params_var('vm');

            %% voltage magnitude violations
            v_Vmin = vm_min - abs(v_);
            v_Vmax = abs(v_) - vm_max;

            %% assemble event function value
            efv = [v_Vmin; v_Vmax];
        end

        function [nx, cx, s] = cpf_callback_vlim(obj, k, nx, cx, px, s, opt, mm, dm, mpopt)
            %% initialize
            if k == 0   %% check for base case voltage violations
                %% convert cx.x back to v_
                ad = mm.get_userdata('aux_data');

                %% get current node voltage magnitudes and bounds
                [v_, ~] = obj.cpf_convert_x(cx.x, ad, 1);
                [~, vm_min, vm_max] = obj.params_var('vm');

                %% violated voltage magnitudes
                if any(abs(v_) < vm_min) || any(abs(v_) > vm_max)
                    %% find node(s) with violated lim(s)
                    ib = find([abs(v_) < vm_min; abs(v_) > vm_max]);
                    nb = length(vm_min);
                    msg = '';
                    for j = 1:length(ib)
                        b = ib(j);          %% index of critical node event of interest
                        if b > nb
                            b = b - nb;
                            nlabel = obj.set_type_label('node', b, dm);
                            msg = sprintf('%snode voltage magnitude limit violated in base case: %s exceeds Vmax limit %g p.u.',...
                               msg, nlabel, vm_max(b));
                        else
                            nlabel = obj.set_type_label('node', b, dm);
                            msg = sprintf('%snode voltage magnitude limit violated in base case: %s exceeds Vmin limit %g p.u.',...
                               msg, nlabel, vm_min(b));
                        end
                    end

                    %% prepare to terminate
                    s.done = 1;
                    s.done_msg = msg;
                end
            end

            %% skip if finalize or done
            if k < 0 || s.done
                return;
            end

            %% handle event
            ev = pne_detected_event(s.events, 'VLIM', 1);   %% zero only
            if ~isempty(ev)
                if opt.verbose > 3
                    msg = sprintf('%s\n    ', ev.msg);
                else
                    msg = '';
                end

                %% find the bus(es) and which lim(s)
                ib = ev.idx;            %% event function index
                [~, vm_min, vm_max] = obj.params_var('vm');
                nb = length(vm_min);
                for j = 1:length(ib)
                    b = ib(j);          %% index of critical node event of interest
                    if b > nb
                        b = b - nb;
                        nlabel = obj.set_type_label('node', b, dm);
                        msg = sprintf('%snode voltage magnitude limit reached\n%s at Vmax limit %g p.u. @ lambda = %.4g, in %d continuation steps',...
                            msg, nlabel, vm_max(b), nx.x(end), k);
                    else
                        nlabel = obj.set_type_label('node', b, dm);
                        msg = sprintf('%snode voltage magnitude limit reached\n%s at Vmin limit %g p.u. @ lambda = %.4g, in %d continuation steps',...
                            msg, nlabel, vm_min(b), nx.x(end), k);
                    end
                end

                %% prepare to terminate
                s.done = 1;
                s.done_msg = msg;
            end
        end

        %%-----  OPF methods  -----
        function [vx_, z_, x_] = opf_convert_x(obj, mmx, ad)
            %% convert (real) math model x to (complex) network model x_
            nv_ = obj.nv / 2;       %% number of voltage vars (sysx=1)
            nz_ = obj.nz;           %% number of state vars
            va = mmx(1:nv_, :);     b = nv_;
            vm = mmx(b+1:b+nv_, :); b = b + nv_;
            zr = mmx(b+1:b+nz_, :); b = b + nz_;
            zi = mmx(b+1:b+nz_, :);
            vx_ = vm .* exp(1j*va);
            z_  = zr+1j*zi;
            if nargout < 2
                vx_ = [vx_; z_];
            elseif nargout > 2
                x_ = [vx_; z_];
            end
        end

        function names = opf_legacy_user_var_names(obj)
            names = {'Va', 'Vm', 'Pg', 'Qg'};
        end
    end     %% methods
end         %% classdef
