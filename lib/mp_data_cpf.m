classdef mp_data_cpf < mp_data
%MP_DATA_CPF  Base class for MATPOWER data model

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
        function obj = mp_data_cpf()
            %% call parent constructor
            obj@mp_data();
            obj.element_classes = ...
                { @dme_bus, @dme_gen, @dme_load_cpf, ...
                    @dme_branch, @dme_shunt_cpf };
        end
    end     %% methods
end         %% classdef
