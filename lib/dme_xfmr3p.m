classdef dme_xfmr3p < dm_element
%DME_XFMR3P  MATPOWER data model class for 3-phase transformer data

%   MATPOWER
%   Copyright (c) 2021-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        fbus    %% bus index vector for "from" port (port 1) (all branches)
        tbus    %% bus index vector for "to" port (port 2) (all branches)
        r       %% series resistance (p.u.) for branches that are on
        x       %% series reactance (p.u.) for branches that are on
        base_kva
        base_kv
%         g_to    %% shunt conductance (p.u.) at "to" end for branches that are on
%         b_fr    %% shunt susceptance (p.u.) at "from" end for branches that are on
%         b_to    %% shunt susceptance (p.u.) at "to" end for branches that are on
%         tm      %% transformer off-nominal turns ratio for branches that are on
%         ta      %% xformer phase-shift angle (radians) for branches that are on
%         rate_a  %% long term flow limit (p.u.) for branches that are on
    end     %% properties

    methods
        function name = name(obj)
            name = 'xfmr3p';
        end

        function label = label(obj)
            label = '3-ph Transformer';
        end

        function label = labels(obj)
            label = '3-ph Transformers';
        end

        function name = cxn_type(obj)
            name = 'bus3p';
        end

        function name = cxn_idx_prop(obj)
            name = {'fbus', 'tbus'};
        end

        function var_names = table_var_names(obj)
            var_names = horzcat( table_var_names@dm_element(obj), ...
                {'bus_fr', 'bus_to', 'r', 'x', 'base_kva', 'base_kv', ...
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
            obj.r = obj.tab.r(obj.on);
            obj.x = obj.tab.x(obj.on);
            obj.base_kva = obj.tab.base_kva(obj.on);
            obj.base_kv  = obj.tab.base_kv(obj.on);
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

        function obj = pp_title(obj, dm, section, out_e, mpopt, fd, varargin)
            if ~strcmp(section, 'det') || strcmp(varargin{1}, 'c_fr')
                %% call parent
                obj = pp_title@dm_element(obj, dm, section, out_e, mpopt, fd, varargin{:});
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
            fprintf(fd, '  %-29s %12.1f kW %12.1f kVAr\n', 'Total 3-ph transformer loss', ...
                sum(sum(ploss(obj.on, :))), sum(sum(qloss(obj.on, :))) );
        end

        function TorF = pp_have_section_det(obj, mpopt, varargin)
            TorF = true;
        end

        function h = pp_get_headers_det(obj, dm, out_e, mpopt, varargin)
            switch varargin{1}
                case 'c_fr'
                    h0 = {'-->  Current Injections at "From" Bus' };
                case 'c_to'
                    h0 = {'', ...
                        '<--  Current Injections at "To" Bus' };
                case 's_fr'
                    h0 = {'', ...
                        '-->  Power Injections at "From" Bus' };
                case 's_to'
                    h0 = {'', ...
                        '<--  Power Injections at "To" Bus' };
            end
            switch varargin{1}(1)
                case 'c'
                    h = {   h0{:}, ...
                            '  3-ph    3-ph Bus  3-ph Bus          Phase A Current  Phase B Current  Phase C Current', ...
                            'Line ID    From ID   To ID    Status   (A)    (deg)     (A)    (deg)     (A)    (deg)', ...
                            '--------  --------  --------  ------  ------  ------   ------  ------   ------  ------' };
                    %%       1234567 123456789 123456789 -----1 123456.89 12345.7 12345.78 12345.7 12345.78 12345.7
                case 's'
                    h = {   h0{:}, ...
                            '  3-ph    3-ph Bus  3-ph Bus          Phase A Power    Phase B Power    Phase C Power', ...
                            'Line ID    From ID   To ID    Status   (kW)   (kVAr)    (kW)   (kVAr)    (kW)   (kVAr)', ...
                            '--------  --------  --------  ------  ------  ------   ------  ------   ------  ------' };
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
