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

        function obj = ext2int(obj, mpopt)
            if ~isfield(obj.mpc, 'order') || obj.mpc.order.state == 'e'
                obj.mpc = ext2int(obj.mpc, mpopt);
            else
%                 warning('mp_data_mpc2/ext2int: data model already in internal format');
            end
        end

        function obj = int2ext(obj, mpopt)
            if isfield(obj.mpc, 'order') && obj.mpc.order.state == 'i'
                obj.mpc = int2ext(obj.mpc, mpopt);
            else
%                 warning('mp_data_mpc2/int2ext: data model already in external format');
            end
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
