classdef mp_branch < mp_element

%   MATPOWER
%   Copyright (c) 2019, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%         name = 'branch';
%     end
    
    methods
        %% constructor
        function obj = mp_branch(varargin)
            obj@mp_element(varargin{:});
            obj.name = 'branch';
            obj.mpc_field = 'branch';
            obj.np = 2;             %% this is a 2 port element
        end

        function obj = build_params(obj, asm, mpc)
%             define_constants;
            branch = mpc.branch;
            nl = obj.nk;

%             f = branch(:, F_BUS);           %% list of "from" buses
%             t = branch(:, T_BUS);           %% list of "to" buses
%             Cf = sparse(f, 1:nl, 1, nb, nl);%% connection matrix for line & from buses
%             Ct = sparse(t, 1:nl, 1, nb, nl);%% connection matrix for line & to buses
% 
%             obj.C = { Cf, Ct };
        end
    end     %% methods
end         %% classdef
