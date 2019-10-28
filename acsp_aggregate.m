classdef acsp_aggregate < ac_aggregate & acsp_model

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        va = [];
        vm = [];
    end
    
    methods
        %% constructor
        function obj = acsp_aggregate(varargin)
            obj@ac_aggregate(varargin{:});
            obj.element_classes = ...
                { @acsp_bus, @acsp_gen, @acsp_load, @acsp_branch };
            if isempty(obj.node)    %% skip if constructed from existing object
                obj.init_set_types();   %% should be called in mp_modeler
                                        %% constructor, if not for:
                                        %% https://savannah.gnu.org/bugs/?52614
            end
        end

        function obj = def_set_types(obj)
            def_set_types@ac_aggregate(obj);        %% call parent first
            obj.set_types.va = 'voltage angle variable';
            obj.set_types.vm = 'voltage magnitude variable';
        end
    end     %% methods
end         %% classdef
