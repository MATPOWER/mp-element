classdef mpe_container < handle
%MPE_CONTAINER  Mix-in class for methods shared by MPE_NETWORK and MP_DATA

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        element_classes %% cell array of function handles of
                        %% constructors for classes for individual
                        %% element types, filled by subclass constructor
        elm_list        %% cell array of dm_element objects
        elm_map         %% key = element name, val = index into elm_list
    end     %% properties

    methods
        function obj = modify_element_classes(obj, class_list)
            %% each element in class_list is either:
            %%  1 - a handle to a constructor to be appended to
            %%      obj.element_classes, or
            %%  2 - a 2-element cell array {A,B} where A is a handle to
            %%      a constructor to replace any element E in the list for
            %%      which isa(E(), B) is true, i.e. B is a char array
            if ~iscell(class_list)
                class_list = {class_list};
            end
            ec = obj.element_classes;   %% list to be updated
            ec0 = {ec{:}};              %% unmodified copy of original list
            for k = 1:length(class_list)
                c = class_list{k};
                if iscell(c)        %% it's a 2-d cell array
                    i = find(cellfun(@(e)isa(e(), c{2}), ec0)); %% find c{2}
                    if ~isempty(i)
                        ec{i} = c{1};                   %% replace with c{1}
                    end
                else                %% it's a single function handle
                    ec{end+1} = c;  %%      append it
                end
            end
            obj.element_classes = ec;
        end

        function elm = elm_by_name(obj, name)
            if isfield(obj.elm_map, name)
                elm = obj.elm_list{obj.elm_map.(name)};
            else
                elm = [];
            end
        end
    end     %% methods
end         %% classdef
