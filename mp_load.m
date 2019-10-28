classdef mp_load < mp_element

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'load';
%     end
    
    methods
        %% constructor
        function obj = mp_load(varargin)
            obj@mp_element(varargin{:});
            obj.name = 'load';
            obj.mpc_field = 'bus';
            obj.np = 1;             %% this is a 1 port element
        end

        function obj = build_params(obj, asm, mpc)
            nb = obj.nk;
            obj.C = { speye(nb) };
        end
    end     %% methods
end         %% classdef
