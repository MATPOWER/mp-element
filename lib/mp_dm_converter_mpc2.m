classdef mp_dm_converter_mpc2 < mp_dm_converter
%MP_DM_CONVERTER_MPC2  MATPOWER data model converter for MATPOWER case v2.

%   MATPOWER
%   Copyright (c) 2021, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

%     properties
%     end     %% properties

    methods
        %% constructor
        function obj = mp_dm_converter_mpc2()
            %% call parent constructor
            obj@mp_dm_converter();
            obj.element_classes = ...
                { @dmce_bus_mpc2, @dmce_gen_mpc2, @dmce_load_mpc2, ...
                    @dmce_branch_mpc2, @dmce_shunt_mpc2 };
        end

        function dm = import(obj, dm, d)
            if ~isstruct(d)
                d = loadcase(d);
            end
            dm = import@mp_dm_converter(obj, dm, d);
        end
    end     %% methods
end         %% classdef
