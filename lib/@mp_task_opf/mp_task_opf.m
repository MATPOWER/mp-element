classdef mp_task_opf < mp_task
%MP_TASK_OPF  MATPOWER task for power flow (OPF).
%   MP_TASK_OPF provides implementation for optimal power flow problem.
%
%   Properties
%       ?
%
%   Methods
%       ?

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        dc      %% true if DC network model (cached in run_pre(), from mpopt)
    end

    methods
        %%-----  constructor  -----
        function obj = mp_task_opf()
            %% call parent constructor
            obj@mp_task();

            obj.tag = 'OPF';
            obj.name = 'Optimal Power Flow';
        end

        %%-----  task methods  -----
        function [m, mpopt] = run_pre(obj, m, mpopt)
            [m, mpopt] = run_pre@mp_task(obj, m, mpopt);     %% call parent

            %% cache DC model flag
            obj.dc = strcmp(upper(mpopt.model), 'DC');

            if ~obj.dc
                alg = upper(mpopt.opf.ac.solver);

                %% check for unsupported solver selection
                if strcmp(alg, 'MINOPF') || strcmp(alg, 'PDIPM') || ...
                        strcmp(alg, 'TRALM') || strcmp(alg, 'SDPOPF')
                    error('mp_task_opf/run_pre: Option ''opf.solver.ac''=''%s'' not supported.', alg);
                end
            end
        end

        %%-----  data model methods  -----
        function dm = data_model_build_post(obj, dm, mpopt)
            dm = data_model_build_post@mp_task(obj, dm, mpopt); %% call parent

            %% pre-process inputs for legacy user vars, constraints, costs
            dm.legacy_user_mod_inputs(mpopt, obj.dc);

            if ~obj.dc
                %% if requested, adjust bus voltage magnitude
                %% limits based on generator Vg setpoint
                use_vg = mpopt.opf.use_vg;
                if use_vg
                    dm.set_bus_v_lims_via_vg(use_vg);
                end
            end
        end

        function dm = data_model_update(obj, mm, nm, dm, mpopt)
            nm.opf_data_model_update(mm, nm, dm, mpopt);
        end

        %%-----  network model methods  -----
        function nm_class = network_model_class_default(obj, dm, mpopt)
            if obj.dc
                nm_class = @mp_network_dc;
            else
                if mpopt.opf.v_cartesian
                    if mpopt.opf.current_balance
                        nm_class = @mp_network_acci;
                    else
                        nm_class = @mp_network_accs;
                    end
                else
                    if mpopt.opf.current_balance
                        nm_class = @mp_network_acpi;
                    else
                        nm_class = @mp_network_acps;
%                        nm_class = @mp_network_acps_test_nln;
                    end
                end
            end
        end

        function nm = network_model_convert_x(obj, mm, nm)
            %% convert solved state from math model to network model soln
            [nm.soln.v, nm.soln.z, nm.soln.x] = nm.opf_convert_x(mm.soln.x, ...
                                                mm.get_userdata('aux_data'));
        end

        %%-----  mathematical model methods  -----
        function mm = math_model_create(obj, nm, dm, mpopt)
            %% switch back to simple opt_model, if possible
            %% I believe opf_model is required for callback functions
            %% that extract mpc from the mm
            mm = opf_model(dm.mpc).init_set_types();
            obj.mm = mm;
        end

        function obj = math_model_add_vars(obj, mm, nm, dm, mpopt)
            nm.opf_add_vars(mm, nm, dm, mpopt);
        end

        function obj = math_model_add_constraints(obj, mm, nm, dm, mpopt)
            nm.opf_add_constraints(mm, nm, dm, mpopt);
        end

        function obj = math_model_add_costs(obj, mm, nm, dm, mpopt)
            nm.opf_add_costs(mm, nm, dm, mpopt);
        end

        function opt = math_model_opt(obj, mm, nm, dm, mpopt)
            if strcmp(mm.problem_type(), 'NLP')
                opt = mpopt2nlpopt(mpopt, mm.problem_type());
            else
                opt = mpopt2qpopt(mpopt, mm.problem_type());
            end

            if obj.dc && strcmp(opt.alg, 'OSQP')
                opt.x0 = [];    %% disable provided starting point for OSQP
            end

            if mpopt.opf.start < 2 && (~obj.dc || ...
                    strcmp(opt.alg, 'MIPS') || strcmp(opt.alg, 'IPOPT'))
                [x0, xmin, xmax] = mm.params_var();     %% init var & bounds
                s = 1;                      %% set init point inside bounds by s
                lb = xmin; ub = xmax;
                lb(xmin == -Inf) = -1e10;   %% replace Inf with numerical proxies
                ub(xmax ==  Inf) =  1e10;   %% temporarily to avoid errors in next line
                x0 = (lb + ub) / 2;         %% set x0 mid-way between bounds
                k = find(xmin == -Inf & xmax < Inf);    %% if only bounded above
                x0(k) = xmax(k) - s;                    %% set just below upper bound
                k = find(xmin > -Inf & xmax == Inf);    %% if only bounded below
                x0(k) = xmin(k) + s;                    %% set just above lower bound

                vv = mm.get_idx();
                bus_dme = dm.elm_by_name('bus');
                gen_nme = nm.elm_by_name('gen');
                if obj.dc
                    Varefs = bus_dme.Va0(find(bus_dme.isref));
                    x0(vv.i1.Va:vv.iN.Va) = Varefs(1);  %% angles set to first reference angle
                else
                    Varefs = bus_dme.Va0(find(bus_dme.isref));
                    Vmax = min(bus_dme.Vmax, 1.5);
                    Vmin = max(bus_dme.Vmin, 0.5);
                    Vm = (Vmax + Vmin) / 2;
                    if mpopt.opf.v_cartesian
                        V = Vm * exp(1j*Varefs(1));
                        x0(vv.i1.Vr:vv.iN.Vr) = real(V);
                        x0(vv.i1.Vi:vv.iN.Vi) = imag(V);
                    else
                        x0(vv.i1.Va:vv.iN.Va) = Varefs(1);  %% angles set to first reference angle
                        x0(vv.i1.Vm:vv.iN.Vm) = Vm;         %% voltage magnitudes
                    end
                end
                if gen_nme.cost.pwl.n > 0
                    gen_dme = dm.elm_by_name('gen');
                    ipwl = gen_nme.cost.pwl.i;
                    maxgc = gen_dme.max_pwl_gencost(ipwl, dm);
                    x0(vv.i1.y:vv.iN.y) = maxgc + 0.1 * abs(maxgc);
                end

                opt.x0 = x0;
            end

            obj.mm_opt = opt;
        end

        function mm = math_model_build_post(obj, mm, nm, dm, mpopt)
            mm = dm.run_userfcn(mm, mpopt);
        end
    end     %% methods
end         %% classdef
