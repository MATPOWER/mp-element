function obj = t_mp_table(quiet)
%T_MP_TABLE  Tests for MP_TABLE (and TABLE).

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargin < 1
    quiet = 0;
end

% define_constants;
if quiet
    verbose = 0;
else
    verbose = 1;
end

nt = 127;
skip_tests_for_tablicious = 1;
table_classes = {@mp_table};
class_names = {'mp_table'};
if have_feature('table')
    table_classes = {@table, table_classes{:}};
    class_names = {'table', class_names{:}};
end
nc = length(table_classes);

t_begin(nc*nt, quiet);

%% set up example data
var1 = [1:6]';
var2 = var1 + var1/10;
var3 = {'one'; 'two'; 'three'; 'four'; 'five'; 'six'};
var4 = 1./var1;
var5 = [3;-5;0;1;-2;4] <= 0;
var_names = {'igr', 'flt', 'str', 'dbl', 'boo'};
var_values = {var1, var2, var3, var4, var5};
row_names = {'row1'; 'second row'; 'row 3'; 'fourth row'; '5'; '6th row'};
dim_names = {'records', 'fields'};

%% run tests
for k = 1:nc
    table_class = table_classes{k};
    cls = upper(class_names{k});

    t = sprintf('%s : constructor - no args : ', cls);
    T = table_class();
    skip = have_feature('octave') && isa(T, 'table') && ...
            skip_tests_for_tablicious;
    t_ok(isa(T, class_names{k}), [t 'class']);
    t_ok(isempty(T), [t 'isempty']);

    t = sprintf('%s : constructor - ind vars : ', cls);
    T = table_class(var1, var2, var3, var4, var5(:));
    t_ok(isa(T, class_names{k}), [t 'class']);
    t_ok(~isempty(T), [t 'not isempty']);
    t_is(size(T), [6 5], 12, [t 'sz = size(T)']);
    t_is(size(T, 1), 6, 12, [t 'sz = size(T, 1)']);
    t_is(size(T, 2), 5, 12, [t 'sz = size(T, 2)']);
    [nr, nz] = size(T);
    t_is([nr, nz], [6 5], 12, [t '[nr, nz] = size(T)']);
    t_ok(isequal(T.Properties.VariableNames, ...
        {'var1', 'var2', 'var3', 'var4', 'Var5'}), [t 'VariableNames'] );
    t_ok(isempty(T.Properties.RowNames), [t 'RowNames'] );
    t_ok(isequal(T.Properties.DimensionNames, {'Row', 'Variables'}), ...
        [t 'DimensionNames'] );
% show_me(T)

    %% subsref () (w/o RowNames)
    t = sprintf('%s : subsref () w/o RowNames : ', cls);
    t_is(size(T([2;4],[1;3;4])), [2 3], 12, [t 'T(ii,jj) : size(T(ii,jj))']);
% show_me(T([2;4],[1;3;4]));
    t_is(size(T(:,1:2:5)), [6 3], 12, [t 'T(:,j1:jN) : size']);
% show_me(T(:,1:2:5));
    t_is(size(T(1:2:5,:)), [3 5], 12, [t 'T(i1:iN,:) : size']);
% show_me(T(1:2:5,:));
    t_is(size(T(5,4)), [1 1], 12, [t 'T(i,j) : size']);
% show_me(T(5,4));

    t = sprintf('%s : constructor - ind vars, w/names : ', cls);
    if skip
        T = table_class(var1, var2, var3, var4, var5, ...
            'VariableNames', var_names, 'RowNames', row_names);
    else
        T = table_class(var1, var2, var3, var4, var5, ...
            'VariableNames', var_names, 'RowNames', row_names, ...
            'DimensionNames', dim_names);
    end
    t_ok(isa(T, class_names{k}), [t 'class'])
    t_ok(~isempty(T), [t 'not isempty']);
    t_is(size(T), [6 5], 12, [t 'sz = size(T)']);
    t_is(size(T, 1), 6, 12, [t 'sz = size(T, 1)']);
    t_is(size(T, 2), 5, 12, [t 'sz = size(T, 2)']);
    [nr, nz] = size(T);
    t_is([nr, nz], [6 5], 12, [t '[nr, nz] = size(T)']);
    t_ok(isequal(T.Properties.VariableNames, var_names), [t 'VariableNames'] );
    t_ok(isequal(T.Properties.RowNames, row_names), [t 'RowNames'] );
    if skip
        t_skip(1, [t 'DimensionNames not yet supported'] );
    else
        t_ok(isequal(T.Properties.DimensionNames, dim_names), [t 'DimensionNames'] );
    end

    %% . indexing
    %% get full variables
    t = sprintf('%s : subsref . : ', cls);
    t_is(T.igr, var1, 12, [t 'T.igr']);
    t_is(T.flt, var2, 12, [t 'T.flt']);
    t_ok(isequal(T.str, var3), [t 'T.str']);
    t_is(T.dbl, var4, 12, [t 'T.dbl']);
    t_is(T.boo, var5, 12, [t 'T.boo']);

    %% get indexed variables
    t_is(T.igr(2), var1(2), 12, [t 'T.igr(i)']);
    t_is(T.flt(4:6), var2(4:6), 12, [t 'T.flt(i1:iN)']);
    t_ok(isequal(T.str{3}, var3{3}), [t 'T.str{i}']);
    t_ok(isequal(T.str(6:-1:4), var3(6:-1:4)), [t 'T.str(iN:-1:i1)']);
    t_is(T.dbl([5;3]), var4([5;3]), 12, [t 'T.dbl(ii)']);
    t_is(T.boo(var5 == 1), var5(var5 == 1), 12, [t 'T.boo(<logical>)']);

    %% set full variables
    t = sprintf('%s : subsasgn . : ', cls);
    T.igr = var1(end:-1:1);
    T.flt = var2(end:-1:1);
    T.str = var3(end:-1:1);
    T.dbl = var4(end:-1:1);
    T.boo = var5(end:-1:1);
    t_is(T.igr, var1(end:-1:1), 12, [t 'T.igr']);
    t_is(T.flt, var2(end:-1:1), 12, [t 'T.flt']);
    t_ok(isequal(T.str, var3(end:-1:1)), [t 'T.str']);
    t_is(T.dbl, var4(end:-1:1), 12, [t 'T.dbl']);
    t_is(T.boo, var5(end:-1:1), 12, [t 'T.boo']);

    %% set indexed variables
    T.igr(2) = 55;
    T.flt(4:6) = [0.4 0.5 0.6];
    T.str{3} = 'tres';
    T.str(6:-1:4) = {'seis'; 'cinco'; 'cuatro'};
    T.dbl([5;3]) = [pi; exp(1)];
    T.boo(var5 == 1) = false;

    v1 = var1(end:-1:1); v1(2) = 55;
    v2 = var2(end:-1:1); v2(4:6) = [0.4 0.5 0.6];
    v3 = var3(end:-1:1);
    v3(6:-1:3) = {'seis'; 'cinco'; 'cuatro'; 'tres'};
    v4 = var4(end:-1:1); v4([5;3]) = [pi; exp(1)];
    v5 = var5(end:-1:1); v5(var5 == 1) = false;
    t_is(T.igr, v1, 12, [t 'T.igr(i)']);
    t_is(T.flt, v2, 12, [t 'T.flt(i1:iN)']);
    t_ok(isequal(T.str, v3), [t 'T.str{i}']);
    t_is(T.dbl, v4, 12, [t 'T.dbl(ii)']);
    t_is(T.boo, v5, 12, [t 'T.boo(<logical>)']);
% show_me(T);

    %% {} indexing
    t = sprintf('%s : subsref {} : ', cls);
    if skip
        t_skip(5, [t 'T{:, j} syntax not yet supported'])
    else
        t_ok(isequal(T{:, 1}, v1), [t 'T{:, 1} == v1']);
        t_ok(isequal(T{:, 2}, v2), [t 'T{:, 2} == v2']);
        t_ok(isequal(T{:, 3}, v3), [t 'T{:, 3} == v3']);
        t_ok(isequal(T{:, 4}, v4), [t 'T{:, 4} == v4']);
        t_ok(isequal(T{:, 5}, v5), [t 'T{:, 5} == v5']);
    end

    t_ok(isequal(T{2, 1}, v1(2)), [t 'T{i, j} == v<j>(i)']);
    t_ok(isequal(T{4, 3}, v3(4)), [t 'T{i, j} == v<j>(i)']);
    %% remove skip when https://github.com/apjanke/octave-tablicious/pull/89
    %% is merged
    if skip
        t_skip(1, [t 'T{end, end} syntax not yet supported'])
    else
        t_ok(isequal(T{end, end}, v5(end)), [t 'T{end, end} == v<end>(end)']);
    end
    if skip
        t_skip(4, [t 'T{ii, j} syntax not yet supported'])
    else
        t_ok(isequal(T{1:3, 2}, v2(1:3)), [t 'T{i1:iN, j} == v<j>(i1:iN)']);
        t_ok(isequal(T{[6;3], 5}, v5([6;3])), [t 'T{ii, j} == v<j>(ii)']);
        t_ok(isequal(T{6:-1:3, [2;4;5]}, [v2(6:-1:3) v4(6:-1:3) v5(6:-1:3)]), [t 'T{iN:-1:i1, jj}']);
        t_ok(isequal(T{:, [4;1]}, [v4 v1]), [t 'T{:, jj} == [v<j1> v<j2> ...]']);
    end

    t = sprintf('%s : subsasgn {} : ', cls);
    T{:, 1} = var1;
    t_is(T.igr, var1, 12, [t 'T{:, 1} = var1']);
    T{:, 2} = var2;
    t_is(T.flt, var2, 12, [t 'T{:, 2} = var2']);
    T{:, 3} = var3;
    t_ok(isequal(T.str, var3), [t 'T{:, 3} = var3']);
    T{:, 4} = var4;
    t_is(T.dbl, var4, 12, [t 'T{:, 4} = var4']);
    T{:, 5} = var5;
    t_is(T.boo, var5, 12, [t 'T{:, 5} = var5']);

    T{2, 1} = v1(2);
    t_ok(isequal(T{2, 1}, v1(2)), [t 'T{i, j} = v<j>(i)']);
    T{4, 3} = v3(4);
    t_ok(isequal(T{4, 3}, v3(4)), [t 'T{i, j} = v<j>(i)']);
    T{1:3, 2} = v2(1:3);
    t_ok(isequal(T.flt(1:3), v2(1:3)), [t 'T{i1:iN, j} = v<j>(i1:iN)']);
    T{[6;3], 5} = v5([6;3]);
    t_ok(isequal(T.boo([6;3]), v5([6;3])), [t 'T{ii, j} = v<j>(ii)']);
    if skip
        t_skip(2, [t 'T{ii, jj} syntax not yet supported'])
    else
        T{6:-1:3, [2;4;5]} = [v2(6:-1:3) v4(6:-1:3) v5(6:-1:3)];
        t_ok(isequal(T{6:-1:3, [2;4;5]}, [v2(6:-1:3) v4(6:-1:3) v5(6:-1:3)]), [t 'T{iN:-1:i1, jj}']);
        T{:, [4;1]} = [v4 v1];
        t_ok(isequal(T{:, [4;1]}, [v4 v1]), [t 'T{:, jj} = [v<j1> v<j2> ...]']);
    end
    T{:, 1} = v1;
    T{:, 2} = v2;
    T{:, 3} = v3;
    T{:, 4} = v4;
    T{:, 5} = v5;

    %% () indexing
    t = sprintf('%s : subsref () : T(ii,jj) : ', cls);
    ii = [2;4]; jj = [1;3;4];
    T2 = T(ii,jj);
% show_me(T2);
    t_is(size(T2), [2 3], 12, [t 'size']);
    t_ok(isequal(T2.Properties.VariableNames, var_names(jj)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v1(ii), v3(ii), v4(ii)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(ii)), [t 'RowNames']);
    if skip
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end

    t = sprintf('%s : subsref () : T(:,j1:s:jN) : ', cls);
    jj = [1:2:5];
    T2 = T(:,1:2:5);
% show_me(T2);
    t_is(size(T2), [6 3], 12, [t 'size']);
    t_ok(isequal(T2.Properties.VariableNames, var_names(jj)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v1, v3, v5}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names), [t 'RowNames']);
    if skip
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end

    t = sprintf('%s : subsref () : T(i1:s:iN,:)) : ', cls);
    ii = [1:2:5];
    T2 = T(1:2:5,:);
% show_me(T2);
    t_is(size(T2), [3 5], 12, [t 'size']);
    t_ok(isequal(T2.Properties.VariableNames, var_names), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v1(ii), v2(ii), v3(ii), v4(ii), v5(ii)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(ii)), [t 'RowNames']);
    if skip
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end

    t = sprintf('%s : subsref () : T(i,j) : ', cls);
    T2 = T(5,4);
% show_me(T2);
    t_is(size(T2), [1 1], 12, [t 'size']);
    t_ok(isequal(T2.Properties.VariableNames, var_names(4)), [t 'VariableNames']);
    t_ok(isequal(table_values(T2), {v4(5)}), [t 'VariableValues']);
    t_ok(isequal(T2.Properties.RowNames, row_names(5)), [t 'RowNames']);
    if skip
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end

    %% remove skip when https://github.com/apjanke/octave-tablicious/pull/89
    %% is merged
    if skip
        t_skip(5, sprintf('%s : T(end, end) not yet supported.', cls));
    else
        t = sprintf('%s : subsref () : T(end, end) : ', cls);
        T2 = T(end, end);
    % show_me(T2);
        t_is(size(T2), [1 1], 12, [t 'size']);
        t_ok(isequal(T2.Properties.VariableNames, var_names(5)), [t 'VariableNames']);
        t_ok(isequal(table_values(T2), {v5(end)}), [t 'VariableValues']);
        t_ok(isequal(T2.Properties.RowNames, row_names(end)), [t 'RowNames']);
        if skip
            t_skip(1, [t 'DimensionNames not yet supported'])
        else
            t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
        end
    end

    if skip
        t_skip(16, sprintf('%s : subsasgn () not yet supported.', cls));
    else
        t = sprintf('%s : subsasgn () : T(ii,jj) : ', cls);
        ii0 = [2;4]; jj = [1;3;4];
        T3 = T(ii0,jj);
        ii = [4;2];
        T(ii,jj) = T3;
        T2 = T(ii,jj);
        t_ok(isequal(T2.Properties.VariableNames, var_names(jj)), [t 'VariableNames']);
        t_ok(isequal(table_values(T2), {v1(ii0), v3(ii0), v4(ii0)}), [t 'VariableValues']);
        t_ok(isequal(T2.Properties.RowNames, row_names(ii)), [t 'RowNames']);
        if skip
            t_skip(1, [t 'DimensionNames not yet supported'])
        else
            t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
        end
        T(ii0, jj) = T2;        %% restore

        t = sprintf('%s : subsasgn () : T(:,j1:s:jN) : ', cls);
        jj = [1:2:5];
        T3 = T(end:-1:1,1:2:5);
        T(:,1:2:5) = T3;
        T2 = T(:,1:2:5);
        t_ok(isequal(T2.Properties.VariableNames, var_names(jj)), [t 'VariableNames']);
        t_ok(isequal(table_values(T2), {v1(end:-1:1), v3(end:-1:1), v5(end:-1:1)}), [t 'VariableValues']);
        t_ok(isequal(T2.Properties.RowNames, row_names), [t 'RowNames']);
        if skip
            t_skip(1, [t 'DimensionNames not yet supported'])
        else
            t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
        end
        T(end:-1:1,1:2:5) = T3;     %% restore

        t = sprintf('%s : subsasgn () : T(i1:s:iN,:)) : ', cls);
        ii0 = [1:2:5];
        T3 = T(1:2:5,:);
        ii = [5:-2:1];
        T(5:-2:1, :) = T3;
        T2 = T(5:-2:1, :);
        t_ok(isequal(T2.Properties.VariableNames, var_names), [t 'VariableNames']);
        t_ok(isequal(table_values(T2), {v1(ii0), v2(ii0), v3(ii0), v4(ii0), v5(ii0)}), [t 'VariableValues']);
        t_ok(isequal(T2.Properties.RowNames, row_names(ii)), [t 'RowNames']);
        if skip
            t_skip(1, [t 'DimensionNames not yet supported'])
        else
            t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
        end
        T(1:2:5, :) = T3;   %% restore

        t = sprintf('%s : subsref () : T(i,j) : ', cls);
        T3 = T(5,4);
        T3{1,1} = exp(1);
        T(5,4) = T3;
        T2 = T(5,4);
        t_ok(isequal(T2.Properties.VariableNames, var_names(4)), [t 'VariableNames']);
        t_ok(isequal(table_values(T2), {exp(1)}), [t 'VariableValues']);
        t_ok(isequal(T2.Properties.RowNames, row_names(5)), [t 'RowNames']);
        if skip
            t_skip(1, [t 'DimensionNames not yet supported'])
        else
            t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
        end
        T{5,4} = pi;
    end

    %% vertical concatenation
    t = sprintf('%s : vertical concatenation : ', cls);
    T1 = T(1:3, :);
    T2 = T(4:6, :);
    T3 = T([4:6 1:3], :);
    T4 = [T2; T1];
    t_ok(isequal(T4, T3), [t '[T1;T2]']);
    T5 = table_class([7;8],[7.7;8.8],{'seven';'eight'},1./[7;8],[-1;2]<=0, ...
        'VariableNames', var_names);
    if skip
        t_skip(1, [t 'RowNames auto-generation not yet supported']);
    else
        T6 = [T5; T2; T1];
        t_ok(isequal(T6, [T5;T3]), [t '[T1;T2;T3]']);
    end

    %% horizontal concatenation
    t = sprintf('%s : horizontal concatenation : ', cls);
    T1 = T(:, 1:3);
    T2 = T(:, 4:5);
    T3 = T(:, [4:5 1:3]);
    T4 = [T2 T1];
    t_ok(isequal(T4, T3), [t '[T1 T2]']);
    T5 = table_class(var3, var2);
    T6 = [T2 T1 T5];
    t_ok(isequal(T6, [T4 T5]), [t '[T1 T2 T3]']);

    %% more {} indexing
    t = sprintf('%s : more subsref {} : ', cls);
    T2 = T(:, [1;2;4;5]);
    T3 = T6(:, 5:6);
    ii = [5;3;1];
    jj = [5:6];
    if skip
        t_skip(6, [t 'T{ii, :} syntax not yet supported'])
    else
        ex = horzcat(v3, var3);
        t_is(T2{ii, :}, [v1(ii) v2(ii) v4(ii) v5(ii)], 12, [t 'T{ii,:} (double)']);
        t_is(T2{:, :}, [v1 v2 v4 v5], 12, [t 'T{:,:} (cell)']);
        t_ok(isequal(T6{ii, 5:6}, ex(ii, :)), [t 'T{ii,j1:j2} (cell)']);
        t_ok(isequal(T6{:, 5:6}, ex), [t 'T{:,j1:j2} (cell)']);
        t_ok(isequal(T3{ii, :}, ex(ii, :)), [t 'T{ii,:} (cell)']);
        t_ok(isequal(T3{:, :}, ex), [t 'T{:,:} (cell)']);
    end

    t = sprintf('%s : more subsasgn {} : ', cls);
    if skip
        t_skip(6, [t 'T{ii, :} syntax not yet supported'])
    else
        T2{ii, :} = [var1(ii) var2(ii) var4(ii) var5(ii)];
        t_is(T2{ii, :}, [var1(ii) var2(ii) var4(ii) var5(ii)], 12, [t 'T{ii,:} (double)']);
        T2{:, :} = [v1 v2 v4 v5];
        t_is(T2{:, :}, [v1 v2 v4 v5], 12, [t 'T{:,:} (cell)']);
        v1 = horzcat(var3, v3);
        v2 = horzcat(v3, var3);
        T6{ii, 5:6} = v1(ii, :);
        t_ok(isequal(T6{ii, 5:6}, v1(ii, :)), [t 'T{ii,j1:j2} (cell)']);
        T6{:, 5:6} = v2;
        t_ok(isequal(T6{:, 5:6}, v2), [t 'T{:,j1:j2} (cell)']);
        T3{ii, :} = v1(ii, :);
        t_ok(isequal(T3{ii, :}, v1(ii, :)), [t 'T{ii,:} (cell)']);
        T3{:, :} = v2;
        t_ok(isequal(T3{:, :}, v2), [t 'T{:,:} (cell)']);
    end

    %% value class, not handle class
    t = sprintf('%s : value class, not handle', cls);
    T2 = T;
    T.igr(2) = 2;
    t_is(T.igr(2), 2, 12, t);
    t_is(T2.igr(2), 55, 12, t);
end

if nargout
    obj = T;
end

t_end;

function Tv = table_values(T)
if isa(T, 'table') && have_feature('matlab')
    Tv = cellfun(@(c)T{:, c}, num2cell(1:size(T, 2)), 'UniformOutput', false);
else
    Tv = T.Properties.VariableValues;
end

function show_me(T)
if isa(T, 'table') && have_feature('octave')
    prettyprint(T)
else
    T
end
