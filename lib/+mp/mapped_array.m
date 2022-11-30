classdef mapped_array < handle
%MP.MAPPED_ARRAY  Cell array indexed by name and numeric index.
%   MA = MP.MAPPED_ARRAY(VALS)
%   MA = MP.MAPPED_ARRAY(VALS, NAMES)
%
%   Input:
%       VALS - cell array of values to be stored
%       NAMES - cell array of names for each element in VALS,
%           where a valid name is any valid variable name, except
%           'p_', 'add_name', 'add_elements', 'delete_elements',
%           'is_index_name', 'name2idx', 'copy' and 'display'
%
%   Currently, arrays are only 1-D.
%
%   Other syntax:
%       val = ma{idx};
%       ma{idx} = val;
%       ma.add_name(idx, name)
%       val = ma.<name>;
%       ma.<name> = val;    %% only if name is already added

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties (Access=private)
        p_ = struct( 'names', {{}}, 'vals', {{}}, 'map', struct() );
    end     %% properties

    methods
        function obj = mapped_array(varargin)
            if nargin
                obj.add_elements(varargin{:});
            end
        end

        function new_obj = copy(obj)
            new_obj = eval(class(obj));  %% create new object

            %% define my_ismethod()
            if have_feature('octave') && have_feature('octave', 'vnum') < 6
                %% Octave 5 ismethod() does not find methods defined in
                %% classdef file, so we use our own here
                my_ismethod = @(o, n)octave5_ismethod_override(o, n);
            else
                my_ismethod = @ismethod;
            end

            %% copy private p_ property which should contain everything
            %% Note: a subclass with additional properties would have to
            %% override copy() if they need special handling
            new_obj.p_ = obj.p_;

            %% make copies of values
            for k = 1:length(obj.p_.vals)
                if isobject(obj.p_.vals{k}) && my_ismethod(obj.p_.vals{k}, 'copy')
                    new_obj.p_.vals{k} = copy(obj.p_.vals{k});
                end
            end
        end

        function len = length(obj)
            len = length(obj.p_.vals);
        end

        function varargout = size(obj, dim)
            varargout = cell(1, nargout);
            w = length(obj.p_.vals);
            h = 1;
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

        function obj = add_names(obj, i0, names)
            if ~iscell(names)
                names = {names};
            end
            N = length(names);
            for k = 1:N
                if isfield(obj.p_.map, names{k}) && obj.p_.map.(names{k}) ~= k
                    error('mp.mapped_array: ''%s'' already refers to element %d', ...
                        names{k}, obj.p_.map.(names{k}));
                end
            end
            for k = N:-1:1
                j = i0 + k - 1;
                obj.p_.names{j} = names{k};
                obj.p_.map.(names{k}) = j;
            end
        end

        function obj = add_elements(obj, vals, names)
            if ~iscell(vals)
                vals = {vals};
            end
            N0 = length(obj.p_.vals);
            N = length(vals);
            if N0
                obj.p_.vals(end+1:end+N) = vals;
                obj.p_.names(end+1:end+N)  = cell(1, N);
            else
                obj.p_.vals = vals;
                obj.p_.names  = cell(1, N);
            end
            if nargin > 2
                add_names(obj, N0+1, names);
            end
        end

        function obj = delete_elements(obj, refs)
            if ~isempty(refs)
                if ischar(refs)
                    refs = {refs};
                end
                if iscell(refs)     %% names
                    name = refs;
                    idx = cellfun(@(x)obj.p_.map.(x), name);
                else                %% idx's
                    idx = refs;
                    name = obj.p_.names(idx);
                end
                [idx, i] = sort(idx, 'descend');
                name = name(i);
                N = length(obj.p_.vals);
                for k = 1:length(idx)
                    %% update map
                    for j = idx(k)+1:N
                        obj.p_.map.(obj.p_.names{j}) = ...
                            obj.p_.map.(obj.p_.names{j}) - 1;
                    end
                    obj.p_.vals(idx(k)) = [];
                    obj.p_.names( idx(k)) = [];
                    obj.p_.map = rmfield(obj.p_.map, name{k});
                    N = N - 1;
                end
            end
        end

        function TorF = is_index_name(obj, name)
            TorF = isfield(obj.p_.map, name);
        end

        function idx = name2idx(obj, name)
            idx = obj.p_.map.(name);
        end

%         function name = idx2name(obj, idx)
%             name = ;
%         end

        function varargout = subsref(obj, s)
            R = length(s) > 1;      %% need to recurse
            switch s(1).type
                case '.'
                    if isfield(obj.p_.map, s(1).subs)   %% name index
                        if R
                            [varargout{1:nargout}] = subsref( ...
                                obj.p_.vals{obj.p_.map.(s(1).subs)}, s(2:end));
                        else
                            varargout{1} = obj.p_.vals{obj.p_.map.(s(1).subs)};
                        end
                    else    %% method calls or properties
                        switch s(1).subs
                            case 'is_index_name'
                                varargout{1} = is_index_name(obj, s(2).subs{:});
                            case 'name2idx'
                                varargout{1} = name2idx(obj, s(2).subs{:});
                            case 'add_names'
                                varargout{1} = add_names(obj, s(2).subs{:});
                            case 'add_elements'
                                varargout{1} = add_elements(obj, s(2).subs{:});
                            case 'delete_elements'
                                varargout{1} = delete_elements(obj, s(2).subs{:});
                            case 'copy'
                                varargout{1} = copy(obj, s(2).subs{:});
                            case 'display'
                                display(obj, s(2).subs{:});
                            case 'p_'   %% properties
                                if R
                                    [varargout{1:nargout}] = subsref(obj.p_, ...
                                        s(2:end));
                                else
                                    varargout{1} = obj.p_;
                                end
                            otherwise   %% unknown methods or properties
                                        %% e.g. defined by subclass
                                if R
                                    [varargout{1:nargout}] = ...
                                        subsref(obj.(s(1).subs), s(2:end));
                                else
                                    varargout{1} = obj.(s(1).subs);
                                end
                        end
                    end
                case '()'
                    if R
                        varargout{1} = subsref(obj.p_.vals(s(1).subs{:}), ...
                            s(2:end))
                    else
                        varargout{1} = obj.p_.vals(s(1).subs{:});
                    end
                case '{}'
                    if R
                        [varargout{1:nargout}] = subsref(obj.p_.vals{s(1).subs{:}}, s(2:end));
                    else
                        %% does not work in Octave for obj{<vector>} because
                        %% is nargout == 0 or 1, while in MATLAB it is
                        %% length(<vector>)
                        %% https://savannah.gnu.org/bugs/index.php?60726
                        [varargout{1:nargout}] = ...
                            obj.p_.vals{s(1).subs{:}};
                    end
            end
        end

        function obj = subsasgn(obj, s, b)
            R = length(s) > 1;      %% need to recurse
            switch s(1).type
                case '.'
                    if isfield(obj.p_.map, s(1).subs)
                        if R
                            obj.p_.vals{ ...
                                obj.p_.map.(s(1).subs)  } = ...
                                    subsasgn(obj.p_.vals{ ...
                                        obj.p_.map.(s(1).subs)  }, ...
                                        s(2:end), b);
                        else
                            obj.p_.vals{ ...
                                obj.p_.map.(s(1).subs)  } = b;
                        end
                    else
                        if R
                            obj.(s(1).subs) = ...
                                subsasgn(obj.(s(1).subs), s(2:end), b);
                        else
                            obj.(s(1).subs) = b;
                        end
                    end
                case '()'
                    if R
                        obj.p_.vals(s(1).subs{:}) = subsasgn(obj.p_.vals(s(1).subs{:}), s(2:end), b);
                    else
                        obj.p_.vals(s(1).subs{:}) = b;
                    end
                case '{}'
                    if R
                        obj.p_.vals{s(1).subs{:}} = ...
                            subsasgn(obj.p_.vals{s(1).subs{:}}, ...
                                s(2:end), b);
                    else
                        obj.p_.vals{s(1).subs{:}} = b;
                    end
            end
        end

        function display(obj)
            if have_feature('matlab')
                display@handle(obj);
            else
                disp(obj), disp('')
            end

            %% find max length for names
            name_len = 6;
            for k = 1:length(obj.p_.vals)
                L = length(obj.p_.names{k});
                if L > name_len
                    name_len = L;
                end
            end

            spc1 = repmat(' ', 1, floor((name_len-2)/2));
            spc2 = repmat(' ', 1, ceil((name_len-2)/2));
            fprintf('   idx   %sname%s   value\n', spc1, spc2);
            ln = repmat('-', 1, name_len+2);
            fprintf('  -----  %s  --------------------\n', ln);
            fmt = sprintf('   %%3d    %%-%ds    %%s\n', name_len);
            for k = 1:length(obj.p_.vals)
                fprintf(fmt, k, obj.p_.names{k}, class(obj.p_.vals{k}))
            end
        end
    end     %% methods
end         %% classdef

function TorF = octave5_ismethod_override(obj, name)
    %OCTAVE5_ISMETHOD_OVERRIDE  Replacement for built-in ISMETHOD() for Octave < 6
    %
    %   In GNU Octave versions before 6.0 the built-in ISMETHOD does not find
    %   methods defined in the classdef file. This implements a replacement
    %   that does. Typical usage:
    %
    %       if have_feature('octave', 'vnum') < 6
    %           my_ismethod = @(o, n)octave5_ismethod_override(o, n);
    %       else
    %           my_ismethod = @ismethod;
    %       end
    %
    %   See also ISMETHOD.

    meths = meta.class.fromName(class(obj)).MethodList;
    for k = 1:length(meths)
        if strcmp(meths{k}.Name, name)
            TorF = true;
            return;
        end
    end
    outTorF = false;
end
