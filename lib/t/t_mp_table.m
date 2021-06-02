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

nt = 38;
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
    skip = have_feature('octave') && strcmp(class_names{k}, 'table') && ...
            skip_tests_for_tablicious;
    table_class = table_classes{k};
    t = sprintf('%s : ', upper(class_names{k}));

    T = table_class();
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
        t_skip(1, [t 'DimensionNames in constructor not implemented'] );
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
    t_is(T.igr(2), var1(2), 12, [t 'get T.igr(2)']);
    t_is(T.flt(4:6), var2(4:6), 12, [t 'get T.flt(4:6)']);
    t_ok(isequal(T.str{3}, var3{3}), [t 'get T.str{3}']);
    t_ok(isequal(T.str(6:-1:4), var3(6:-1:4)), [t 'get T.str(6:-1:4)']);
    t_is(T.dbl([5;3]), var4([5;3]), 12, [t 'get T.dbl([5;3])']);
    t_is(T.boo(var5 == 1), var5(var5 == 1), 12, [t 'get T.boo(var5 == 1)']);

    %% set full variables
    if skip
        T.igr(:) = var1(end:-1:1);
        T.flt(:) = var2(end:-1:1);
        T.str(:) = var3(end:-1:1);
        T.dbl(:) = var4(end:-1:1);
        T.boo(:) = var5(end:-1:1);
    else
        T.igr = var1(end:-1:1);
        T.flt = var2(end:-1:1);
        T.str = var3(end:-1:1);
        T.dbl = var4(end:-1:1);
        T.boo = var5(end:-1:1);
    end
    t_is(T.igr, var1(end:-1:1), 12, [t 'set T.igr']);
    t_is(T.flt, var2(end:-1:1), 12, [t 'set T.flt']);
    t_ok(isequal(T.str, var3(end:-1:1)), [t 'set T.str']);
    t_is(T.dbl, var4(end:-1:1), 12, [t 'set T.dbl']);
    t_is(T.boo, var5(end:-1:1), 12, [t 'set T.boo']);

    %% set indexed variables
    T.igr(2) = 55;
    T.flt(4:6) = [0.4 0.5 0.6];
    if ~skip
        T.str{3} = 'tres';
    end
    T.str(6:-1:4) = {'seis'; 'cinco'; 'cuatro'};
    T.dbl([5;3]) = [pi; exp(1)];
    T.boo(var5 == 1) = false;

    v1 = var1(end:-1:1); v1(2) = 55;
    v2 = var2(end:-1:1); v2(4:6) = [0.4 0.5 0.6];
    v3 = var3(end:-1:1);
    if ~skip, v3{3} = 'tres'; end
    v3(6:-1:4) = {'seis'; 'cinco'; 'cuatro'};
    v4 = var4(end:-1:1); v4([5;3]) = [pi; exp(1)];
    v5 = var5(end:-1:1); v5(var5 == 1) = false;
    t_is(T.igr, v1, 12, [t 'set T.igr(2)']);
    t_is(T.flt, v2, 12, [t 'set T.flt(4:6)']);
    t_ok(isequal(T.str, v3), [t 'set T.str{3}']);
    t_is(T.dbl, v4, 12, [t 'set T.dbl([5;3])']);
    t_is(T.boo, v5, 12, [t 'set T.boo(var5 == 1)']);
% show_me(T)
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
