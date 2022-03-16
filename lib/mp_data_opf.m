classdef mp_data_opf < mp_data
%MP_DATA_OPF  Base class for MATPOWER data model

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
    end     %% properties

    methods
        %% constructor
        function obj = mp_data_opf()
            %% call parent constructor
            obj@mp_data();
            obj.element_classes = ...
                { @dme_bus_opf, @dme_gen_opf, @dme_load_opf, ...
                    @dme_branch_opf, @dme_shunt_opf, ...
                    @dme_bus3p_opf, @dme_gen3p_opf, @dme_load3p_opf, ...
                    @dme_line3p_opf, @dme_xfmr3p_opf, ...
                    @dme_buslink_opf };
        end

        function [out, add] = pp_flags(obj, mpopt)
            %% call parent
            [out, add_] = pp_flags@mp_data(obj, mpopt);
            suppress = add_.suppress;
            s0 = add_.s0;

            %% add limit options
            out_s = struct( ...
                'all', 0, ...
                'any', 0, ...
                'elm', s0 );

            if out.all == 1
                out_s.all = 2;
            elseif out.all == -1 && ~suppress
                out_s.all = mpopt.out.lim.all;
            end

            if out_s.all == -1
                out_s.elm.bus     = ~suppress && mpopt.out.lim.v;
                out_s.elm.branch  = ~suppress && mpopt.out.lim.line;
                out_s.elm.gen     = ~suppress && ...
                                    (mpopt.out.lim.pg || mpopt.out.lim.qg);
            else
                out_s.elm.bus     = out_s.all;
                out_s.elm.branch  = out_s.all;
                out_s.elm.gen     = out_s.all;
            end
            out_s.any = any(cell2mat( ...
                    cellfun(@double, struct2cell(out_s.elm), ...
                            'UniformOutput', false) ));

            out.sec.lim = out_s;
            out.any = out.any || out_s.any;

            %% return additional data
            if nargout > 1
                add = add_;
            end
        end

        function sections = pp_section_list(obj, out)
            sections = pp_section_list@mp_data(obj, out);
            sections{end+1} = 'lim';
        end

        function h = pp_get_headers_other(obj, section, out_s, mpopt)
            switch section
                case 'lim'
                    h = {};
                otherwise
                    error('mp_data_opf:pp_get_headers_other: unknown section ''%s''', section);
            end
        end
    end     %% methods
end         %% classdef
