classdef mp_data_mpc2 < mp_data
    properties
        mpc
    end     %% properties

    methods
        %% constructor
        function obj = mp_data_mpc2(m)
            %% call parent constructor
            obj@mp_data();
            obj.mpc = loadcase(m);
        end

        function display(obj)
            class(obj)
            obj.mpc
        end

        function print_soln(obj, fname)
        end

        function save_soln(obj, fname)
        end
    end     %% methods
end         %% classdef
