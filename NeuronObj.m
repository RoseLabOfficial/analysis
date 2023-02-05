classdef NeuronObj
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        membrane_capcitance
        input_resistance
        resting_potential
        exctiatory_reversal_potential
        inhibitory_reversal_potential
        current_clamps
        stimulus
    end

    methods
        function obj = NeuronObj(datafile, paradigms, stats_duration)
            if nargin > 0
                [m, n] = size(stats_duration);
                obj(n, m) = obj;
            end
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end