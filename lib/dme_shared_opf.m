classdef dme_shared_opf < handle
%DME_SHARED_OPF  Mix-in class for MATPOWER data model elements for MP_DATA_OPF

%   MATPOWER
%   Copyright (c) 2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        ctol    %% constraint violation tolerance
        ptol    %% shadow price tolerance
    end     %% properties

    methods
        function obj = pp_set_tols_lim(obj, mpopt)
            obj.ctol = mpopt.opf.violation;
            obj.ptol = 1e-4;
        end

        function TorF = pp_have_section_other(obj, section, mpopt, varargin)
            switch section
                case 'lim'
                    TorF = obj.pp_have_section_lim(section, mpopt, varargin{:});
                otherwise
                    error('dme_shared_opf:pp_have_section_other: unknown section ''%s''', section);
            end
        end

        function rows = pp_rows_other(obj, dm, section, out_e, mpopt, varargin)
            switch section
                case 'lim'
                    if isempty(obj.ptol)
                        obj.pp_set_tols_lim(mpopt);
                    end
                    rows = obj.pp_rows_lim(dm, out_e, mpopt, varargin{:});
                otherwise
                    error('dme_shared_opf:pp_rows_other: unknown section ''%s''', section);
            end
        end

        function str = pp_title_str_other(obj, section, mpopt, varargin)
            switch section
                case 'lim'
                    str = obj.pp_title_str_lim(mpopt, varargin{:});
                otherwise
                    error('dme_shared_opf:pp_title_str_other: unknown section ''%s''', section);
            end
        end

        function h = pp_get_headers_other(obj, dm, section, out_e, mpopt, varargin)
            switch section
                case 'lim'
                    h = obj.pp_get_headers_lim(dm, out_e, mpopt, varargin{:});
                otherwise
                    error('dme_shared_opf:pp_get_headers_other: unknown section ''%s''', section);
            end
        end

        function obj = pp_data_other(obj, dm, section, rows, out_e, mpopt, fd, varargin)
            switch section
                case 'lim'
                    obj.pp_data_lim(dm, rows, out_e, mpopt, fd, varargin{:});
                otherwise
                    error('dme_shared_opf:pp_data_other: unknown section ''%s''', section);
            end
        end

        function TorF = pp_have_section_lim(obj, mpopt, varargin)
            TorF = false;
        end

        function rows = pp_rows_lim(obj, dm, out_e, mpopt, varargin)
            if out_e == 2       %% all rows
                rows = -1;
            elseif out_e == 1   %% binding rows
                rows = obj.pp_binding_rows_lim(dm, out_e, mpopt, varargin{:});
            else                %% no rows
                rows = 0;
            end
        end

        function rows = pp_binding_rows_lim(obj, dm, out_e, mpopt, varargin)
            rows = 0;           %% no rows
        end

        function str = pp_title_str_lim(obj, mpopt, varargin)
            str = sprintf('%s Constraints', obj.label);
        end

        function h = pp_get_headers_lim(obj, dm, out_e, mpopt, varargin)
            h = {};
        end

        function obj = pp_data_lim(obj, dm, rows, out_e, mpopt, fd, varargin)
            if ~isempty(rows) && rows(1) == -1  %% all rows
                for k = 1:obj.nr
                    fprintf(fd, '%s\n', ...
                        obj.pp_data_row_lim(dm, k, out_e, mpopt, fd, varargin{:}));
                end
            else
                for k = 1:length(rows)
                    fprintf(fd, '%s\n', ...
                        obj.pp_data_row_lim(dm, rows(k), out_e, mpopt, fd, varargin{:}));
                end
            end
        end

        function str = pp_data_row_lim(obj, dm, k, out_e, mpopt, fd, varargin)
            str = '';
        end
    end     %% methods
end         %% classdef
