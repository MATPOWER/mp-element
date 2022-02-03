classdef mp_math_cpf_acp < mp_math_cpf
%MP_MATH_CPF_ACP  MATPOWER mathematical model for continuation power flow (CPF) problem.
%   ?
%
%   MP_MATH_CPF_ACP ... power flow ...
%
%   Properties
%       ? - ?
%
%   Methods
%       ?

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
        %% constructor
        function obj = mp_math_cpf_acp()
            obj@mp_math_cpf();
            obj.element_classes = { mme_bus_pf_acp, mme_gen_pf_ac, ...
                mme_branch_pf_ac, ...
                @mme_bus3p, @mme_gen3p, @mme_line3p, @mme_xfmr3p, ...
                @mme_buslink_pf_acp };
        end

        function efv = event_vlim(obj, cx, opt, nm, dm, mpopt)
            %% convert cx.x back to v_
            %% get current node voltage magnitudes and bounds
            [v_, ~] = obj.cpf_convert_x(cx.x, nm, 1);
            [~, vm_min, vm_max] = nm.params_var('vm');

            %% voltage magnitude violations
            v_Vmin = vm_min - abs(v_);
            v_Vmax = abs(v_) - vm_max;

            %% assemble event function value
            efv = [v_Vmin; v_Vmax];
        end

        function [nx, cx, s] = callback_vlim(obj, k, nx, cx, px, s, opt, nm, dm, mpopt)
            %% initialize
            if k == 0   %% check for base case voltage violations
                %% convert cx.x back to v_
                %% get current node voltage magnitudes and bounds
                [v_, ~] = obj.cpf_convert_x(cx.x, nm, 1);
                [~, vm_min, vm_max] = nm.params_var('vm');

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
                            nlabel = nm.set_type_label('node', b, dm);
                            msg = sprintf('%snode voltage magnitude limit violated in base case: %s exceeds vm upper bound %g p.u.',...
                               msg, nlabel, vm_max(b));
                        else
                            nlabel = nm.set_type_label('node', b, dm);
                            msg = sprintf('%snode voltage magnitude limit violated in base case: %s exceeds vm lower bound %g p.u.',...
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
                [~, vm_min, vm_max] = nm.params_var('vm');
                nb = length(vm_min);
                for j = 1:length(ib)
                    b = ib(j);          %% index of critical node event of interest
                    if b > nb
                        b = b - nb;
                        nlabel = nm.set_type_label('node', b, dm);
                        msg = sprintf('%snode voltage magnitude limit reached\n%s at vm upper bound %g p.u. @ lambda = %.4g, in %d continuation steps',...
                            msg, nlabel, vm_max(b), nx.x(end), k);
                    else
                        nlabel = nm.set_type_label('node', b, dm);
                        msg = sprintf('%snode voltage magnitude limit reached\n%s at vm lower bound %g p.u. @ lambda = %.4g, in %d continuation steps',...
                            msg, nlabel, vm_min(b), nx.x(end), k);
                    end
                end

                %% prepare to terminate
                s.done = 1;
                s.done_msg = msg;
            end
        end
    end     %% methods
end         %% classdef
