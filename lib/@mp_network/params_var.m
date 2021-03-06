function [v0, vl, vu, vt] = params_var(obj, vtype, name, idx)
%PARAMS_VAR  Returns initial value, lower bound and upper bound for opt variables.
%   [V0, VL, VU] = OBJ.PARAMS_VAR(VTYPE)
%   [V0, VL, VU] = OBJ.PARAMS_VAR(VTYPE, NAME)
%   [V0, VL, VU] = OBJ.PARAMS_VAR(VTYPE, NAME, IDX_LIST)
%   [V0, VL, VU, VT] = OBJ.PARAMS_VAR(...)
%   Returns the initial value V0, lower bound VL and upper bound VU for
%   the full optimization variable vector, or for a specific named or named
%   and indexed variable set. Optionally also returns a corresponding char
%   vector VT of variable types, where 'C', 'I' and 'B' represent continuous
%   integer and binary variables, respectively.
%
%   Examples:
%       [x0, xmin, xmax] = obj.params_var();
%       [Pg, Pmin, Pmax] = obj.params_var('Pg');
%       [zij0, zijmin, zijmax, ztype] = obj.params_var('z', {i, j});
%   
%   See also OPT_MODEL/PARAMS_VAR.

%   MATPOWER
%   Copyright (c) 2008-2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

if nargout > 3
    have_vt = 1;
else
    have_vt = 0;
end
var = obj.(vtype);
if nargin < 3
    v0 = []; vl = []; vu = []; vt = char([]);
    %% calls to substruct() are relatively expensive, so we pre-build the
    %% structs for addressing cell and numeric array fields, updating only
    %% the subscripts before use
    sc = struct('type', {'.', '{}'}, 'subs', {'', 1});  %% cell array field
    for k = 1:var.NS
        name = var.order(k).name;
        idx = var.order(k).idx;
        if isempty(idx)
            v0 = [ v0; var.data.v0.(name) ];
            vl = [ vl; var.data.vl.(name) ];
            vu = [ vu; var.data.vu.(name) ];
            if have_vt
                N = var.idx.N.(name);
                vt0 = var.data.vt.(name);
                if isscalar(vt0) && N > 1 
                    vt = [ vt char(vt0 * ones(1, N)) ];
                else
                    vt = [ vt vt0 ];
                end
            end
        else
            % (calls to substruct() are relatively expensive ...
            % sc = substruct('.', name, '{}', idx);
            % ... so replace it with these more efficient lines)
            sc(1).subs = name;
            sc(2).subs = idx;
            v0 = [ v0; subsref(var.data.v0, sc) ];
            vl = [ vl; subsref(var.data.vl, sc) ];
            vu = [ vu; subsref(var.data.vu, sc) ];
            if have_vt
                % (calls to substruct() are relatively expensive ...
                % sn = substruct('.', name, '()', idx);
                % ... so replace it with these more efficient lines)
                sn = sc; sn(2).type = '()';
                N = subsref(var.idx.N, sn);
                vt0 = subsref(var.data.vt, sc);
                if isscalar(vt0) && N > 1 
                    vt = [ vt char(vt0 * ones(1, N)) ];
                else
                    if ~isempty(vt0)
                        vt = [ vt vt0 ];
                    end
                end
            end
        end
    end
else
    if isfield(var.idx.N, name)
        if nargin < 4 || isempty(idx)
            v0 = var.data.v0.(name);
            vl = var.data.vl.(name);
            vu = var.data.vu.(name);
            if have_vt
                N = var.idx.N.(name);
                vt0 = var.data.vt.(name);
                if isscalar(vt0) && N > 1 
                    vt = char(vt0 * ones(1, N));
                else
                    vt = vt0;
                end
            end
        else
            % (calls to substruct() are relatively expensive ...
            % sc = substruct('.', name, '{}', idx);
            % ... so replace it with these more efficient lines)
            sc = struct('type', {'.', '{}'}, 'subs', {name, idx});
            v0 = subsref(var.data.v0, sc);
            vl = subsref(var.data.vl, sc);
            vu = subsref(var.data.vu, sc);
            if have_vt
                % (calls to substruct() are relatively expensive ...
                % sn = substruct('.', name, '()', idx);
                % ... so replace it with these more efficient lines)
                sn = sc; sn(2).type = '()';
                N = subsref(var.idx.N, sn);
                vt0 = subsref(var.data.vt, sc);
                if isscalar(vt0) && N > 1 
                    vt = char(vt0 * ones(1, N));
                else
                    vt = vt0;
                end
            end
        end
    else
        v0 = [];
        vl = [];
        vu = [];
        if have_vt
            vt = [];
        end
    end
end
