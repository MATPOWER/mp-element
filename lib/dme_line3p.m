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
            label = '3-ph Line';
        end

        function label = labels(obj)
            label = '3-ph Lines';
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
%             vars = {'pl1_fr', 'ql1_fr', ...
%                     'pl2_fr', 'ql2_fr', ...
%                     'pl3_fr', 'ql3_fr', ...
%                     'pl1_to', 'ql1_to', ...
%                     'pl2_to', 'ql2_to', ...
%                     'pl3_to', 'ql3_to' };
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

        function pretty_print(obj, dm, section, out_e, mpopt, fd, varargin)
            switch section
                case 'det'
                    %% compute currents
                    s_fr = [  obj.tab.pl1_fr + 1j * obj.tab.ql1_fr ...
                            obj.tab.pl2_fr + 1j * obj.tab.ql2_fr ...
                            obj.tab.pl3_fr + 1j * obj.tab.ql3_fr    ];
                    s_to = [  obj.tab.pl1_to + 1j * obj.tab.ql1_to ...
                            obj.tab.pl2_to + 1j * obj.tab.ql2_to ...
                            obj.tab.pl3_to + 1j * obj.tab.ql3_to    ];
                    t = dm.elements.bus3p.tab;    %% bus3p table
                    vm = [ t.vm1 t.vm2 t.vm3 ] .* (t.base_kv/sqrt(3) * [1 1 1]);
                    va = [ t.va1 t.va2 t.va3 ];
                    v_ = vm .* exp(1j * va * pi/180);
                    i_fr = conj( s_fr ./ v_(obj.fbus, :));
                    i_to = conj( s_to ./ v_(obj.tbus, :));
                    c_fr = struct('cm', abs(i_fr), 'ca', angle(i_fr) * 180/pi);
                    c_to = struct('cm', abs(i_to), 'ca', angle(i_to) * 180/pi);

                    pretty_print@dm_element(obj, dm, section, out_e, mpopt, fd, 'c_fr', c_fr, varargin{:});
                    pretty_print@dm_element(obj, dm, section, out_e, mpopt, fd, 'c_to', c_to, varargin{:});
                    pretty_print@dm_element(obj, dm, section, out_e, mpopt, fd, 's_fr', s_fr, varargin{:});
                    pretty_print@dm_element(obj, dm, section, out_e, mpopt, fd, 's_to', s_to, varargin{:});
                otherwise
                    %% call parent
                    pretty_print@dm_element(obj, dm, section, out_e, mpopt, fd, varargin{:});
            end
        end

        function TorF = pp_have_section_sum(obj, mpopt, varargin)
            TorF = true;
        end

        function obj = pp_data_sum(obj, dm, rows, out_e, mpopt, fd, varargin)
            %% call parent
            pp_data_sum@dm_element(obj, dm, rows, out_e, mpopt, fd, varargin{:});

            %% print generation summary
            t = obj.tab;
            ploss = [ t.pl1_fr t.pl2_fr t.pl3_fr ] + ...
                    [ t.pl1_to t.pl2_to t.pl3_to ];
            qloss = [ t.ql1_fr t.ql2_fr t.ql3_fr ] + ...
                    [ t.ql1_to t.ql2_to t.ql3_to ];
            fprintf(fd, '  %-29s %12.1f kW %12.1f kVAr\n', 'Total 3-ph line loss', ...
                sum(sum(ploss(obj.on, :))), sum(sum(qloss(obj.on, :))) );
        end

        function TorF = pp_have_section_det(obj, mpopt, varargin)
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, varargin)
            if strcmp(varargin{1}, 'c_fr')
                h1 = pp_get_headers_det@dm_element(obj, dm, out_e, mpopt, varargin{:});
            else
                h1 = {};
            end
            switch varargin{1}
                case 'c_fr'
                    h2 = {'-->  Current Injections at "From" Bus' };
                case 'c_to'
                    h2 = {'', ...
                        '<--  Current Injections at "To" Bus' };
                case 's_fr'
                    h2 = {'', ...
                        '-->  Power Injections at "From" Bus' };
                case 's_to'
                    h2 = {'', ...
                        '<--  Power Injections at "To" Bus' };
            end
            switch varargin{1}(1)
                case 'c'
                    h = [ h1 h2 ...
                        {   '  3-ph    3-ph Bus  3-ph Bus          Phase A Current  Phase B Current  Phase C Current', ...
                            'Line ID    From ID   To ID    Status   (A)    (deg)     (A)    (deg)     (A)    (deg)', ...
                            '--------  --------  --------  ------  ------  ------   ------  ------   ------  ------' } ];
                    %%       1234567 123456789 123456789 -----1 123456.89 12345.7 12345.78 12345.7 12345.78 12345.7
                case 's'
                    h = [ h1 h2 ...
                        {   '  3-ph    3-ph Bus  3-ph Bus          Phase A Power    Phase B Power    Phase C Power', ...
                            'Line ID    From ID   To ID    Status   (kW)   (kVAr)    (kW)   (kVAr)    (kW)   (kVAr)', ...
                            '--------  --------  --------  ------  ------  ------   ------  ------   ------  ------' } ];
                    %%       1234567 123456789 123456789 -----1 1234567.9 12345.7 123456.8 12345.7 123456.8 12345.7
            end
        end

        function str = pp_data_row_det(obj, dm, k, out_e, mpopt, fd, varargin)
            switch varargin{1}(1)
                case 'c'
                    cm = varargin{2}.cm(k, :);
                    ca = varargin{2}.ca(k, :);
                    str = sprintf('%7d %9d %9d %6d %9.2f %7.1f %8.2f %7.1f %8.2f %7.1f', ...
                        obj.tab.uid(k), obj.tab.bus_fr(k), obj.tab.bus_to(k), ...
                        obj.tab.status(k), ...
                        cm(:, 1), ca(:, 1), ...
                        cm(:, 2), ca(:, 2), ...
                        cm(:, 3), ca(:, 3) );
                case 's'
                    sk = varargin{2}(k, :);
                    str = sprintf('%7d %9d %9d %6d %9.1f %7.1f %8.1f %7.1f %8.1f %7.1f', ...
                        obj.tab.uid(k), obj.tab.bus_fr(k), obj.tab.bus_to(k), ...
                        obj.tab.status(k), ...
                        real(sk(:, 1)), imag(sk(:, 1)), ...
                        real(sk(:, 2)), imag(sk(:, 2)), ...
                        real(sk(:, 3)), imag(sk(:, 3)) );
            end
        end
    end     %% methods
end         %% classdef
