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

nt = 75;
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
    t = sprintf('%s : ', cls);

    T = table_class();
    skip = have_feature('octave') && isa(T, 'table') && ...
            skip_tests_for_tablicious;
    t_ok(isa(T, class_names{k}), [t 'constructor - no args'])

    T = table_class(var1, var2, var3, var4, var5(:));
    t_ok(isa(T, class_names{k}), [t 'constructor - ind vars'])
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

    %% () indexing (no RowNames)
    t_is(size(T([2;4],[1;3;4])), [2 3], 12, [t 'size(T(ii,jj))']);
% show_me(T([2;4],[1;3;4]));
    t_is(size(T(:,1:2:5)), [6 3], 12, [t 'size(T(:,j1:jN))']);
% show_me(T(:,1:2:5));
    t_is(size(T(1:2:5,:)), [3 5], 12, [t 'size(T(i1:iN,:))']);
% show_me(T(1:2:5,:));
    t_is(size(T(5,4)), [1 1], 12, [t 'size(T(i,j))']);
% show_me(T(5,4));

    if skip
        T = table_class(var1, var2, var3, var4, var5, ...
            'VariableNames', var_names, 'RowNames', row_names);
    else
        T = table_class(var1, var2, var3, var4, var5, ...
            'VariableNames', var_names, 'RowNames', row_names, ...
            'DimensionNames', dim_names);
    end
    t_ok(isa(T, class_names{k}), [t 'constructor - ind vars, w/names'])
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

    %% get full variables
    t_is(T.igr, var1, 12, [t 'get T.igr']);
    t_is(T.flt, var2, 12, [t 'get T.flt']);
    t_ok(isequal(T.str, var3), [t 'get T.str']);
    t_is(T.dbl, var4, 12, [t 'get T.dbl']);
    t_is(T.boo, var5, 12, [t 'get T.boo']);

    %% get indexed variables
    t_is(T.igr(2), var1(2), 12, [t 'get T.igr(i)']);
    t_is(T.flt(4:6), var2(4:6), 12, [t 'get T.flt(i1:iN)']);
    t_ok(isequal(T.str{3}, var3{3}), [t 'get T.str{i}']);
    t_ok(isequal(T.str(6:-1:4), var3(6:-1:4)), [t 'get T.str(iN:-1:i1)']);
    t_is(T.dbl([5;3]), var4([5;3]), 12, [t 'get T.dbl(ii)']);
    t_is(T.boo(var5 == 1), var5(var5 == 1), 12, [t 'get T.boo(<logical>)']);

    %% set full variables
    T.igr = var1(end:-1:1);
    T.flt = var2(end:-1:1);
    T.str = var3(end:-1:1);
    T.dbl = var4(end:-1:1);
    T.boo = var5(end:-1:1);
    t_is(T.igr, var1(end:-1:1), 12, [t 'set T.igr']);
    t_is(T.flt, var2(end:-1:1), 12, [t 'set T.flt']);
    t_ok(isequal(T.str, var3(end:-1:1)), [t 'set T.str']);
    t_is(T.dbl, var4(end:-1:1), 12, [t 'set T.dbl']);
    t_is(T.boo, var5(end:-1:1), 12, [t 'set T.boo']);

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
    t_is(T.igr, v1, 12, [t 'set T.igr(i)']);
    t_is(T.flt, v2, 12, [t 'set T.flt(i1:iN)']);
    t_ok(isequal(T.str, v3), [t 'set T.str{i}']);
    t_is(T.dbl, v4, 12, [t 'set T.dbl(ii)']);
    t_is(T.boo, v5, 12, [t 'set T.boo(<logical>)']);
% show_me(T);

    %% {} indexing
    t = sprintf('%s {} indexing : ', cls);
    if skip
        t_skip(5, [t 'T{:, j} syntax not yet supported'])
    else
        t_ok(isequal(T{:, 1}, v1), [t 'T{:, 1} = v1']);
        t_ok(isequal(T{:, 2}, v2), [t 'T{:, 2} = v2']);
        t_ok(isequal(T{:, 3}, v3), [t 'T{:, 3} = v3']);
        t_ok(isequal(T{:, 4}, v4), [t 'T{:, 4} = v4']);
        t_ok(isequal(T{:, 5}, v5), [t 'T{:, 5} = v5']);
    end

    t_ok(isequal(T{2, 1}, v1(2)), [t 'T{i, j} = v<j>(i)']);
    t_ok(isequal(T{4, 3}, v3(4)), [t 'T{i, j} = v<j>(i)']);
    if skip
        t_skip(5, [t 'T{ii, j} syntax not yet supported'])
    else
        t_ok(isequal(T{1:3, 2}, v2(1:3)), [t 'T{i1:iN, j} = v<j>(i1:iN)']);
        t_ok(isequal(T{[6;3], 5}, v5([6;3])), [t 'T{ii, j} = v<j>(ii)']);
        t_ok(isequal(T{:, 4}, v4), [t 'T{:, j} = v<j>']);
        t_ok(isequal(T{6:-1:3, [2;4;5]}, [v2(6:-1:3) v4(6:-1:3) v5(6:-1:3)]), [t 'T{iN:-1:i1, jj}']);
        t_ok(isequal(T{:, [4;1]}, [v4 v1]), [t 'T{:, jj} = [v<j1> v<j2> ...]']);
    end

    %% () indexing
    t = sprintf('%s () indexing : ', cls);
    t = sprintf('%s : T2 = T(ii,jj) : ', cls);
    ii = [2;4]; jj = [1;3;4];
    T2 = T(ii,jj);
    t_is(size(T2), [2 3], 12, [t 'size(T2)']);
    t_ok(isequal(T2.Properties.VariableNames, var_names(jj)), [t 'VariableNames']);
    if skip
        t_skip(1, [t 'T{:, j} syntax not yet supported'])
    else
        t_ok(isequal({T2{:,1}, T2{:,2}, T2{:,3}}, {v1(ii), v3(ii), v4(ii)}), [t 'values']);
    end
    if skip
        t_skip(1, [t 'DimensionNames not yet supported'])
    else
        t_ok(isequal(T2.Properties.DimensionNames, dim_names), [t 'DimensionNames']);
    end
    t_ok(isequal(T2.Properties.RowNames, row_names(ii)), [t 'RowNames']);

    t_is(size(T([2;4],[1;3;4])), [2 3], 12, [t 'size(T(ii,jj))']);
% show_me(T([2;4],[1;3;4]));
    t_is(size(T(:,1:2:5)), [6 3], 12, [t 'size(T(:,j1:s:jN))']);
% show_me(T(:,1:2:5));
    t_is(size(T(1:2:5,:)), [3 5], 12, [t 'size(T(i1:s:iN,:))']);
% show_me(T(1:2:5,:));
    t_is(size(T(5,4)), [1 1], 12, [t 'size(T(i,j))']);
% show_me(T(5,4));

    %% vertical concatenation
    t = 'vertical concatenation : ';
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
    t = 'horizontal concatenation : ';
    T1 = T(:, 1:3);
    T2 = T(:, 4:5);
    T3 = T(:, [4:5 1:3]);
    T4 = [T2 T1];
    t_ok(isequal(T4, T3), [t '[T1 T2]']);
    T5 = table_class(var3, var2);
    T6 = [T2 T1 T5];
    t_ok(isequal(T6, [T4 T5]), [t '[T1 T2 T3]']);

    %% more {} indexing
    t = sprintf('%s more {} indexing : ', cls);
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

    %% value class, not handle class
    t = 'value class, not handle';
    T2 = T;
    T.igr(2) = 2;
    t_is(T.igr(2), 2, 12, t);
    t_is(T2.igr(2), 55, 12, t);
end

if nargout
    obj = T;
end

t_end;

function show_me(T)
if isa(T, 'table') && have_feature('octave')
    prettyprint(T)
else
    T
end
