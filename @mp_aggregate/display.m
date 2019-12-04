function display(obj)
%DISPLAY  Displays the object.
%   Called when semicolon is omitted at the command-line. Displays the details
%   of the variables, constraints, costs included in the model.
%
%   See also MP_MODELER.

%   MATPOWER
%   Copyright (c) 2008-2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%% base element info
display@mp_element(obj)
fprintf('\n');

%% nodes and states
obj.display_set('node', 'NODES');
obj.display_set('state', 'STATES');

%% variables
vvars = obj.model_vvars();
zvars = obj.model_zvars();
for k = 1:length(vvars)
    obj.display_set(vvars{k}, obj.set_types.(vvars{k}));
end
for k = 1:length(zvars)
    obj.display_set(zvars{k}, obj.set_types.(zvars{k}));
end

%% elements
model_params = obj.model_params();
fprintf('ELEMENTS\n')
fprintf('========\n')
fprintf('       name          N      np    nz    class, param(m,n))\n');
fprintf('  ------------   --------  ----  ----  --------------------\n');
for k = 1:length(obj.mpe_list)
    mpe = obj.mpe_list{k};
    fprintf(' %11s %11d %5d %5d    %s', mpe.name, mpe.nk, mpe.np, mpe.nz, class(mpe));
    
    for j = 1:length(model_params)
        pn = model_params{j};   %% parameter name
        if ~isempty(mpe.(pn))
            [m, n] = size(mpe.(pn));
            fprintf(', %s(%d,%d)', pn, m, n);
        end
    end
    fprintf('\n');
%     mpe
end
