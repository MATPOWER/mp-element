classdef dme_branch < dm_element
%DME_BRANCH  Abstract base class for MATPOWER data model branch table

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        fbus    %% bus index vector for "from" port (port 1) (all branches)
        tbus    %% bus index vector for "to" port (port 2) (all branches)
        R       %% series resistance (p.u.) for branches that are on
        X       %% series reactance (p.u.) for branches that are on
        B       %% series reactance (p.u.) for branches that are on
        tap     %% transformer off-nominal turns ratio for branches that are on
        shift   %% xformer phase-shift angle (radians) for branches that are on
        rate_a  %% long term flow limit (p.u.) for branches that are on
    end     %% properties

    methods
        %% constructor
        function obj = dme_branch()
            obj@dm_element();   %% call parent constructor
            obj.name = 'branch';
            obj.table = 'branch';
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus_fr', 'bus_to', 'r', 'x', 'g_fr', 'b_fr', ...
                'g_to', 'b_to', 'sm_ub_a', 'sm_ub_b', 'sm_ub_c', ...
                'cm_ub_a', 'cm_ub_b', 'cm_ub_c', 'vad_lb', 'vad_ub', ...
                'tm', 'ta', ... %% remove these when we separate out xformers
                'pl_fr', 'ql_fr', 'pl_to', 'ql_to', ...
                'mu_flow_fr_ub', 'mu_flow_to_ub', ...
                'mu_vad_lb', 'mu_vad_ub'});
%                 'mu_sm_fr_ub', 'mu_sm_to_ub', ...
%                 'mu_pl_fr_ub', 'mu_pl_to_ub', ...
%                 'mu_cm_fr_ub', 'mu_cm_to_ub', ...
        end

        function obj = initialize(obj, dm)
            initialize@dm_element(obj, dm);     %% call parent

            %% get bus mapping info
            b2i = dm.elements.bus.ID2i;     %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            obj.fbus = b2i(obj.tab.bus_fr);
            obj.tbus = b2i(obj.tab.bus_to);
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elements.bus.status;    %% bus status

            %% update status of branches connected to isolated/offline buses
            obj.status = obj.status & bs(obj.fbus) & ...
                                      bs(obj.tbus);

            %% call parent to fill in on/off
            update_status@dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            obj.R  = obj.tab.r(obj.on);
            obj.X  = obj.tab.x(obj.on);
            %%-----  HACK ALERT  -----
            %% use the from and to properties properly in 
            obj.B  = obj.tab.b_fr(obj.on) + obj.tab.b_to(obj.on);
            %%-----  end of HACK  -----
            obj.tap    = obj.tab.tm(obj.on);
            obj.shift  = obj.tab.ta(obj.on) * pi/180;
            obj.rate_a = obj.tab.sm_ub_a(obj.on) / dm.baseMVA;
        end

        function obj = update(obj, dm, varargin)
            %% obj.update(dm, name1, val1, name2, val2, ...)
            %% obj.update(dm, idx, name1, val1, name2, val2, ...)

            %% define named indices into data matrices
            [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
                TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
                ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
            baseMVA = dm.baseMVA;

            n = length(varargin);
            if rem(n, 2)    %% odd
                idx = obj.on(varargin{1});
                s = 2;      %% starting arg index
            else            %% even
                idx = obj.on;
                s = 1;      %% starting arg index
            end
            for k = s:2:n-1
                val = varargin{k+1};
                switch varargin{k}
                    case 'Sf'
                        obj.tab.pl_fr(idx) = real(val) * baseMVA;
                        obj.tab.ql_fr(idx) = imag(val) * baseMVA;
                    case 'St'
                        obj.tab.pl_to(idx) = real(val) * baseMVA;
                        obj.tab.ql_to(idx) = imag(val) * baseMVA;
                    case 'Pf'
                        obj.tab.pl_fr(idx) = val * baseMVA;
                    case 'Pt'
                        obj.tab.pl_to(idx) = val * baseMVA;
                    case 'Qf'
                        obj.tab.ql_fr(idx) = val * baseMVA;
                    case 'Qt'
                        obj.tab.ql_to(idx) = val * baseMVA;
                    case {'muSf', 'muPf'}
                        obj.tab.mu_flow_fr_ub(idx) = val / baseMVA;
                    case {'muSt', 'muPt'}
                        obj.tab.mu_flow_to_ub(idx) = val / baseMVA;
                    case 'muAngmin'
                        obj.tab.mu_vad_lb(idx) = val * pi/180;
                    case 'muAngmax'
                        obj.tab.mu_vad_ub(idx) = val * pi/180;
                end
            end
        end

        %%-----  OPF methods  -----
        function [A, l, u, i] = opf_branch_ang_diff_params(obj, dm, ignore)
            %% from makeAang()
            nb = dm.elements.bus.n;
            branch = obj.tab;

            if ignore
                A  = sparse(0, nb);
                l  = [];
                u  = [];
                i  = [];
            else
                i = find( ...
                    branch.vad_lb ~= 0 & ...
                        (branch.vad_lb > -360 | branch.vad_ub == 0) | ...
                    branch.vad_ub ~= 0 & ...
                        (branch.vad_ub <  360 | branch.vad_lb == 0) );
                n = length(i);

                if n > 0
                    ii = [(1:n)'; (1:n)'];
                    jj = [obj.fbus(i); obj.tbus(i)];
                    A = sparse(ii, jj, [ones(n, 1); -ones(n, 1)], n, nb);
                    l = branch.vad_lb(i);
                    u = branch.vad_ub(i);
                    l(l < -360) = -Inf;
                    u(u >  360) =  Inf;
                    l = l * pi/180;
                    u = u * pi/180;
                else
                    A = sparse(0, nb);
                    l =[];
                    u =[];
                end
            end
        end
    end     %% methods
end         %% classdef
