classdef mp_model < handle
%MP_MODEL  MATPOWER Model abstract base class.
%   Each concrete MATPOWER Element class must inherit, at least indirectly,
%   from both MP_ELEMENT and MP_MODEL.
%
%   MP_MODEL provides propoerties and methods related to the specific
%   model and formulation (e.g. DC version, AC polar power version, etc.)
%
%   Properties
%       subclasses provide properties for model parameters
%
%   Methods
%       model_name() - returns string w/name of model/formulation
%       model_tag() - returns string w/short label for model/formulation
%       model_params() - cell array of names of model parameters

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         mp_model_field = '';
%     end
    
    methods
        function name = model_name(obj)
            error('model_name() method not implemented');
        end
        function tag = model_tag(obj)
            error('model_tag() method not implemented');
        end
        function params = model_params(obj)
            error('model_params() method not implemented');
        end

        function model_class = find_model_class(obj)
            if isa(obj, 'mp_model')
                tab = obj.superclass_tab({'mp_model', 'mp_element'});

                %% select among classes that are not mp_element's ...
                j = find(tab.ii(:, 2) == 0);
                %% ... the mp_model class with longest inheritance path
                %% back to mp_model
                [~, k] = max(tab.ii(j, 1));

                assert(~isempty(k));
                model_class = tab.name{j(k)};
            else
                model_class = '<none>';
            end
        end

        function [tab, ii] = superclass_tab(obj, roots, mcls, tab, ii, level)
            %% obj - an mp_model object
            %% roots - cell array of names of root classes to search for
            %% mcls - meta class object (undocumented MATLAB and Octave) for
            %%        the class of interest
            %% tab - struct with fields:
            %%      name - cell array of class names
            %%      ii   - matrix where each row corresponds to entry in name
            %%             field, and column j corresponds to j-th entry in
            %%             roots, value is generations of inheritance required
            %%             to reach that root
            %% ii - row of tab.ii corresponding to mcls on input (before
            %%      traversing parent classes)
            %% level - indent level for this class when printing
            %%         (hard-coded 'verbose' variable not 0)
            verbose = 0;
            if nargin < 6
                level = 0;
                if nargin < 5
                    ii = [];
                    if nargin < 4
                        tab = struct('name', {{}}, 'ii', []);
                        if nargin < 3
                            mcls = meta.class.fromName(class(obj));
                            if nargin < 2
                                roots = {};
                            end
                        end
                    end
                end
            end
            n = length(roots);
            if isempty(ii)
                ii = zeros(1, n);
            end
            done = zeros(1, n);
            i0 = ii;
            for k = 1:n
                if strcmp(mcls.Name, roots{k})
                    done(k) = 1;
                end
            end
            if verbose
                prefix = repmat(' ', 1, 4*level);
                fprintf('%s %s\n', prefix, mcls.Name);
            end
            if have_fcn('octave')
                sclist = mcls.SuperClassList;
            else
                sclist = {};
                for k = 1:length(mcls.SuperclassList)
                    sclist{end+1} = mcls.SuperclassList(k);
                end
            end
            for k = 1:length(sclist)
                [tab, iii] = obj.superclass_tab(roots, sclist{k}, tab, i0, level+1);
                for k = 1:n
                    if done(k) == 1
                        ii(k) = 1;
                    elseif iii(k) && i0(k) < iii(k) + 1
                        ii(k) = iii(k) + 1;
                    end
                end
            end
            if verbose
                fmt = repmat(' %d', 1, n);
                fprintf(['%s' fmt '\n'], prefix, ii);
            end
            tab.name{end+1, 1} = mcls.Name;
            tab.ii(end+1, 1:n) = ii;
        end
    end     %% methods
end         %% classdef
