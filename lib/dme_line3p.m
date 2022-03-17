classdef dme_line3p < dm_element
%DME_LINE3P  MATPOWER data model class for 3-phase line data

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        fbus    %% bus index vector for "from" port (port 1) (all lines)
        tbus    %% bus index vector for "to" port (port 2) (all lines)
        freq    %% system frequency, in Hz
        lc      %% index into lc_tab for lines that are on
        len     %% length for lines that are on
        lc_tab  %% line construction table
        ys      %% cell array of 3x3 series admittance matrices for lc rows
        yc      %% cell array of 3x3 shunt admittance matrices for lc rows
    end     %% properties

    methods
        function name = name(obj)
            name = 'line3p';
        end

        function label = label(obj)
            label = 'Line (3-ph)';
        end

        function label = labels(obj)
            label = 'Lines (3-ph)';
        end

        function name = cxn_type(obj)
            name = 'bus3p';
        end

        function name = cxn_idx_prop(obj)
            name = {'fbus', 'tbus'};
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus_fr', 'bus_to', 'lc', 'len', ...
                 'pl1_fr', 'ql1_fr', 'pl2_fr', 'ql2_fr', 'pl3_fr', 'ql3_fr', ...
                 'pl1_to', 'ql1_to', 'pl2_to', 'ql2_to', 'pl3_to', 'ql3_to' ...
                 });
        end

%         function vars = export_vars(obj, task)
%             switch task
%                 case 'PF'
%                     vars = {'ql_to', 'pl_to', 'ql_fr', 'pl_fr'};
%                 case 'CPF'
%                     vars = {'ql_to', 'pl_to', 'ql_fr', 'pl_fr'};
%                 case 'OPF'
%                     vars = {'ql_to', 'pl_to', 'ql_fr', 'pl_fr', ...
%                         'mu_flow_fr_ub', 'mu_flow_to_ub', ...
%                         'mu_vad_lb', 'mu_vad_ub'};
%                 otherwise
%                     vars = 'all';
%             end
%         end

        function obj = initialize(obj, dm)
            initialize@dm_element(obj, dm); %% call parent

            %% get bus mapping info
            b2i = dm.elements.bus3p.ID2i;   %% bus num to idx mapping

            %% set bus index vectors for port connectivity
            obj.fbus = b2i(obj.tab.bus_fr);
            obj.tbus = b2i(obj.tab.bus_to);
        end

        function obj = update_status(obj, dm)
            %% get bus status info
            bs = dm.elements.bus3p.tab.status;  %% bus status

            %% update status of branches connected to isolated/offline buses
            obj.tab.status = obj.tab.status & bs(obj.fbus) & ...
                                              bs(obj.tbus);

            %% call parent to fill in on/off
            update_status@dm_element(obj, dm);
        end

        function obj = build_params(obj, dm)
            nlc = size(obj.lc, 1);
            obj.ys = zeros(nlc, 6);
            obj.yc = zeros(nlc, 6);
            obj.lc = obj.tab.lc(obj.on);
            obj.len = obj.tab.len(obj.on);

            %% build Ys and Yc for relevant lines
            idx = unique(obj.lc);
            rr = obj.lc_tab.r(idx, :);
            xx = obj.lc_tab.x(idx, :);
            cc = obj.lc_tab.c(idx, :);
%             rr = obj.lc_tab{idx, 2:7};
%             xx = obj.lc_tab{idx, 8:13};
%             cc = obj.lc_tab{idx, 14:19};
            for k = 1:length(idx)
                R = obj.vec2symmat(rr(k, :));
                X = obj.vec2symmat(xx(k, :));
                C = obj.vec2symmat(cc(k, :));
                Ys = inv(R + 1j * X);
                Yc = 1j * 2*pi * obj.freq * 1e-9 * C;
                obj.ys(k, :) = obj.symmat2vec(Ys);
                obj.yc(k, :) = obj.symmat2vec(Yc);
            end
        end

        function M = vec2symmat(obj, v)
            % makes a symmetric matrix from a vector of 6 values
            M = [v(1) v(2) v(3);
                 v(2) v(4) v(5);
                 v(3) v(5) v(6) ];
        end

        function v = symmat2vec(obj, M)
            % extracts a vector of 6 values from a matrix assumed to be symmetric
            v = [M(1, :) M(2,2:3) M(3,3)];
        end
    end     %% methods
end         %% classdef
