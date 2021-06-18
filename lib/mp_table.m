classdef mp_table
%MP_TABLE  Very basic TABLE class for Octave or older Matlab.

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        Properties = struct('VariableNames',  {{}}, ...
                            'VariableValues', {{}}, ...
                            'RowNames',       {{}}, ...
                            'DimensionNames', {{}}    );
    end     %% properties

    methods
        function obj = mp_table(varargin)
            args = varargin;

            if nargin
                %% extract named arguments
                [var_names, row_names, dim_names, args] = ...
                    extract_named_args(obj, args);

                %% set default variable names
                nv = length(args);          %% number of variables
                if length(var_names) < nv
                    for k = nv:-1:length(var_names)+1
                        var_names{k} = inputname(k);
                        if isempty(var_names{k})
                            var_names{k} = sprintf('Var%d', k);
                        end
                    end
                end

                %% set default dimension names
                if isempty(dim_names)
                    dim_names = {'Row', 'Variables'};
                end

                %% check dimensions
                if length(var_names) > nv
                    warning('mp_table: ignoring %d extra VariableNames provided', length(var_names)-nv);
                end
                nr = size(args{1}, 1);      %% number of rows
                for k = 2:nv
                    if size(args{k}, 1) ~= nr
                        error('mp_table: variable %d expected to have %d rows', k, nr);
                    end
                end

                %% set table properties
                obj.Properties.VariableNames = var_names;
                obj.Properties.DimensionNames = dim_names;
                obj.Properties.RowNames = row_names;
                obj.Properties.VariableValues = args;
                obj.Properties.Map = cell2struct(num2cell(1:nv), var_names, 2);
            end
        end

        function new_obj = copy(obj)
            new_obj = eval(class(obj));  %% create new object

            %% make copies of properties
            if have_feature('octave')
                s1 = warning('query', 'Octave:classdef-to-struct');
                warning('off', 'Octave:classdef-to-struct');
            end
            props = fieldnames(obj);
            if have_feature('octave')
                warning(s1.state, 'Octave:classdef-to-struct');
            end
            for k = 1:length(props)
                new_obj.(props{k}) = obj.(props{k});
            end

            %% make copies of values
            for k = 1:length(obj.Properties.VariableValues)
                if isobject(obj.Properties.VariableValues{k}) && ...
                        ismethod(obj.Properties.VariableValues{k}, 'copy')
                    new_obj.Properties.VariableValues{k} = ...
                        copy(obj.Properties.VariableValues{k});
                end
            end
        end

        function varargout = size(obj, dim)
            varargout = cell(1, nargout);
            w = length(obj.Properties.VariableValues);
            if w
                h = size(obj.Properties.VariableValues{1}, 1);
            else
                h = 0;
            end
            if nargin == 2
                if dim == 1
                  varargout{1} = h;
                elseif dim == 2
                  varargout{1} = w;
                else
                  varargout{1} = 1;
                end
            elseif nargout == 0 || nargout == 1
                varargout{1} = [h, w];
            else
                varargout{1} = h;
                varargout{2} = w;
                [varargout{3:end}] = deal(1);
            end
        end

        function b = subsref(obj, s)
            switch s(1).type
                case '.'
                    if strcmp(s(1).subs, 'Properties')
                        b = obj.Properties;
                    else
                        p = obj.Properties;
                        b = p.VariableValues{p.Map.(s(1).subs)};
                    end
                case '()'
                    b = obj(s(1).subs{:});
                case '{}'
                    b = obj{s(1).subs{:}};
            end
            if length(s) > 1    %% recurse
                b = subsref(b, s(2:end));
            end
        end

        function obj = subsasgn(obj, s, b)
            R = length(s) > 1;      %% need to recurse
            switch s(1).type
                case '.'
                    if strcmp(s(1).subs, 'Properties')
                        if R
                            obj.Properties = ...
                                subsasgn(obj.Properties, s(2:end), b);
                        else
                            obj.Properties = b;
                        end
                    else
                        idx = obj.Properties.Map.(s(1).subs);
                        if R
                            obj.Properties.VariableValues{idx} = ...
                                subsasgn(obj.Properties.VariableValues{idx}, ...
                                    s(2:end), b);
                        else
                            obj.Properties.VariableValues{idx} = b;
                        end
                    end
                case '()'
                    if R
                        obj(s(1).subs{:}) = subsasgn(obj(s(1).subs{:}), s(2:end), b);
                    else
                        obj(s(1).subs{:}) = b;
                    end
                case '{}'
                    if R
                        obj{s(1).subs{:}} = subsasgn(obj{s(1).subs{:}}, s(2:end), b);
                    else
                        obj{s(1).subs{:}} = b;
                    end
            end
        end

        function [var_names, row_names, dim_names, args] = ...
                extract_named_args(obj, args)
            var_names = {};
            row_names = {};
            dim_names = {};
            for k = length(args)-1:-1:1
                if ischar(args{k})
                    switch args{k}
                        case 'VariableNames'
                            var_names = args{k+1};
                            args(k:k+1) = [];
                        case 'RowNames'
                            row_names = args{k+1};
                            args(k:k+1) = [];
                        case 'DimensionNames'
                            dim_names = args{k+1};
                            args(k:k+1) = [];
                    end
                end
            end
        end

        function display(obj)
            allrows = 1;

            [nr, nv] = size(obj);

            %% rows to display
            max_nr = 25;
            first_last = 10;
            if nr <= max_nr || allrows
                rr = [1:nr]';
            else
                rr = [1:first_last nr-first_last+1:nr]';
            end

            %% var names
            vn = obj.Properties.VariableNames;
            rn = obj.Properties.RowNames;

            %% col_widths
            if isempty(rn)
                cwr = 0;
            else
                cwr = max(cellfun(@length, rn));    %% col width for row names
            end
            for v = nv:-1:1
                cw(v) = length(vn{v});  %% actual column width
                cwv(v) = 0;             %% column width for displayed values
                var = obj.Properties.VariableValues{v};
                if iscell(var)
                    for r = 1:length(rr)
                        s = sprintf('''%s''', var{rr(r)});
                        if length(s) > cw(v);
                            cw(v) = length(s);
                        end
                        if length(s) > cwv(v);
                            cwv(v) = length(s);
                        end
                    end
                else
                    for r = 1:length(rr)
                        s = sprintf('%g', var(rr(r)));
                        if length(s) > cw(v);
                            cw(v) = length(s);
                        end
                        if length(s) > cwv(v);
                            cwv(v) = length(s);
                        end
                    end
                end
                spc{v} = repmat(' ', 1, floor((cw(v)-cwv(v))/2));
            end

            fprintf('\n %dx%d %s\n\n', nr, nv, class(obj));

            %% display var names & lines
            fprintf('  %s', repmat(' ', 1, cwr));
            for v = 1:nv
                fmt = sprintf('  %%-%ds', cw(v));
                name = sprintf('%s', vn{v});
                spcs = repmat(' ', 1, floor((cw(v)-length(name))/2));
                fprintf(fmt, [spcs name]);
            end
            fprintf('\n');
            fprintf('  %s', repmat(' ', 1, cwr));
            for v = 1:nv
                s = repmat('-', 1, cw(v));
                fprintf('  %s', s);
            end
            fprintf('\n');

            %% display var values
            for r = 1:length(rr)
                %% row names
                if ~isempty(rn)
                    fmt = sprintf('  %%-%ds', cwr);
                    fprintf(fmt, rn{rr(r)});
                end

                %% variable values
                for v = 1:nv
                    var = obj.Properties.VariableValues{v};
                    if iscell(var)
                        val = sprintf('%s''%s''', spc{v}, var{rr(r), 1});
                        fmt = sprintf('  %%-%ds', cw(v));
                    else
                        val = sprintf('%g%s', var(rr(r), 1), spc{v});
                        fmt = sprintf('  %%%ds', cw(v));
                    end
                    fprintf(fmt, val);
                end
                fprintf('\n');

                %% print continuation dots if not displaying all rows
                if nr > length(rr) && r == first_last
                    for k = 1:3
                        for v = 1:nv
                            spcs = repmat(' ', 1, floor(cwv(v)/2));
                            if iscell(var)
                                dot = sprintf('%s%s.', spc{v}, spcs);
                                fmt = sprintf('  %%-%ds', cw(v));
                            else
                                dot = sprintf('.%s%s', spcs, spc{v});
                                fmt = sprintf('  %%%ds', cw(v));
                            end
                            fprintf(fmt, dot);
                        end
                        fprintf('\n');
                    end
                end
            end
            fprintf('\n');
        end
    end     %% methods
end         %% classdef
