classdef mp_data_opf < mp_data
%MP_DATA_OPF  Base class for MATPOWER data model

%   MATPOWER
%   Copyright (c) 2020-2022, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
    end     %% properties

    methods
        %% constructor
        function obj = mp_data_opf()
            %% call parent constructor
            obj@mp_data();
            obj.element_classes = ...
                { @dme_bus_opf, @dme_gen_opf, @dme_load, ...
                    @dme_branch_opf, @dme_shunt, ...
                    @dme_bus3p, @dme_gen3p, @dme_load3p, ...
                    @dme_line3p, @dme_xfmr3p, ...
                    @dme_buslink };
        end
    end     %% methods
end         %% classdef
