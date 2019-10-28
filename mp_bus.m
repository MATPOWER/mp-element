classdef mp_bus < mp_element

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'bus';
%     end

    methods
        %% constructor
        function obj = mp_bus(varargin)
            obj@mp_element(varargin{:});
            obj.name = 'bus';
            obj.mpc_field = 'bus';
            obj.np = 1;             %% this is a 1 port element
        end

        function obj = add_nodes(obj, asm, mpc)
            %% define constants
            [PQ, PV, REF, NONE, BUS_I] = idx_bus;
    
            nb = size(mpc.bus, 1);      %% number of buses
            asm.add_node(obj.name, nb, mpc.bus(:, BUS_I));
        end
    end     %% methods
end         %% classdef
