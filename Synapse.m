classdef Synapse

    properties
        conductance
        current
        reversal_potential
    end

    methods
        function obj = Synapse(reversal_potential)
            obj.conductance = zeros(1, 1);
            obj.current = zeros(1, 1);
            obj.reversal_potential = reversal_potential;
        end

        function obj = put_conductance(obj, conductance)
            obj.conductance = conductance;
        end
    
        function obj = call(obj, membrane_potential)
             obj.current = obj.conductance.*(membrane_potential - obj.reversal_potential);
        end
    end

end