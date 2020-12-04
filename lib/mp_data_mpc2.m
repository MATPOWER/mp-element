classdef mp_data_mpc2 < mp_data
%MP_DATA_MPC2  Implementation of MATPOWER data model for MATPOWER case format v2

%   MATPOWER
%   Copyright (c) 2020, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://matpower.org for more info.

    properties
        mpc
    end     %% properties

    methods
        %% constructor
        function obj = mp_data_mpc2(m)
            %% call parent constructor
            obj@mp_data();
            obj.element_classes = ...
                { @dme_bus_mpc2, @dme_gen_mpc2, @dme_branch_mpc2 };

            if nargin
                %% load case and create mappings
                obj.mpc = loadcase(m);
                obj.create_model();
            end
        end

        %%-----  HACK ALERT  -----
        %% This should not be necessary once we get all of the elements
        %% (e.g. loads, shunts, etc.) defined in their own classes
        function n = online(obj, name, tab_name)
            n = online@mp_data(obj, name);      %% call parent
            
            if n == 0 && isfield(obj.mpc, tab_name)
                n = size(obj.mpc.(tab_name), 1);
            end
        end
        %%-----  end of HACK  -----

        function obj = ext2int(obj, mpopt)
            if ~isfield(obj.mpc, 'order') || obj.mpc.order.state == 'e'
                if nargin > 1
                    obj.mpc = ext2int(obj.mpc, mpopt);
                else
                    obj.mpc = ext2int(obj.mpc);
                end
            else
%                 warning('mp_data_mpc2/ext2int: data model already in internal format');
            end
        end

        function obj = int2ext(obj, mpopt)
            if isfield(obj.mpc, 'order') && obj.mpc.order.state == 'i'
                if nargin > 1
                    obj.mpc = int2ext(obj.mpc, mpopt);
                else
                    obj.mpc = int2ext(obj.mpc);
                end
            else
%                 warning('mp_data_mpc2/int2ext: data model already in external format');
            end
        end

        function display(obj)
            fprintf('Data Model class : %s\n', class(obj));
            mpc = obj.mpc
        end

        function print_soln(obj, fname)
        end

        function save_soln(obj, fname)
        end
    end     %% methods
end         %% classdef
