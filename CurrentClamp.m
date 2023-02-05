classdef CurrentClamp

    properties
        injected_current
        activation_potential
        steady_state_potential
        spikes_per_stimulus
        alpha_multiplier
        beta_multiplier
        membrane_potential
        times
        sampling_rate
    end

    properties
        alpha
        beta
        depolarization
        hyperpolarization
        leakage_current
        membrane_current
        active_current
        excitatory_current
        inhibitory_current
    end

    methods
        function obj = CurrentClamp(data, parameters)
            [m, n] = size(parameters.Iinj);
            obj(m, n) = obj;
            for i = 1: 1: m
                for j = 1: 1: n
                    obj(i, j).injected_current = parameters.Iinj(i, j);
                    obj(i, j).activation_potential = parameters.Eact(i, j);
                    obj(i, j).steady_state_potential = parameters.Ess(i, j);
                    obj(i, j).spikes_per_stimulus = parameters.sps(i, j);
                    obj(i, j).alpha_multiplier = parameters.xalpha(i, j);
                    obj(i, j).beta_multiplier = parameters.xbeta(i, j);
                    obj(i, j).membrane_potential = data(:, j+1);
                    obj(i, j).times = data.times;
                    obj(i, j).sampling_rate = 1/(obj(i, j).times(2, 1) - obj(i, j).times(1, 1));
                    obj(i, j).alpha = 0;
                    obj(i, j).beta = 0;
                    obj(i, j).depolarization = zeros(size(obj(i, j).times));
                    obj(i, j).hyperpolarization = zeros(size(obj(i, j).times));
                    obj(i, j).leakage_current = zeros(size(obj(i, j).times));
                    obj(i, j).membrane_current = zeros(size(obj(i, j).times));
                    obj(i, j).active_current = zeros(size(obj(i, j).times));
                    obj(i, j).excitatory_current = zeros(size(obj(i, j).times));
                    obj(i, j).inhibitory_current = zeros(size(obj(i, j).times));
                end
            end
        end

        function obj = compute_leakage_current(obj, input_resistance)
            [m, n] = size(obj);
            for i = 1: 1: m
                for j = 1: 1: n
                    obj(i, j).leakage_current = (1/input_resistance).*(obj(i, j).membrane_potential - obj(i, j).resting_potential);
                end
            end
        end

        function obj = compute_membrane_current(obj, membrane_capacitance)
            [m, n] = size(obj);
            for i = 1: 1: m
                for j = 1: 1: n
                    Vm_appended = cat(1, obj(i, j).membrane_potential(1, :), obj(i, j).membrane_potential, obj(i, j).membrane_potential(end, :));
                    dvdt = diff(Vm_appended);
                    obj(i, j).membrane_current = membrane_capacitance.*dvdt(1:end-1);
                end
            end
        end

        function obj = compute_active_conductances(obj, input_resistance, threshold_potential)
            [m, n] = size(obj);
            for i = 1: 1: m
                for j = 1: 1: n
                    obj(i, j).alpha = (1/input_resistance).*(obj(i, j).activation_potential - obj(i, j).steady_state_potential);
                    obj(i, j).beta = obj(i, j).alpha.*(threshold_potential - obj(i, j).steady_state_potential);
                end
            end
        end

        function obj = compute_active_current(obj, resting_potential, threshold_potential)
            [m, n] = size(obj);
            for i = 1: 1: m
                for j = 1: 1: n
                    obj(i, j).active_current = (obj(i, j).alpha_multiplier.*obj(i, j).alpha.* ...
                        (obj(i, j).membrane_potential - threshold_potential).* ...
                        (obj(i, j).membrane_potential - resting_potential)) + ...
                        (obj(i, j).beta_multiplier.*obj(i, j).beta.*(obj(i, j).membrane_potential - resting_potential));
                end
            end
        end

        function obj = compute_excitatory_current(obj, excitatory_conductance, excitatory_reversal_potential)
            [m, n] = size(obj);
            for i = 1: 1: m
                for j = 1: 1: n
                    obj(i, j).excitatory_current = excitatory_conductance.*(obj(i, j).membrane_potential - excitatory_reversal_potential);
                end
            end
        end

        function obj = compute_inhibitory_current(obj, inhibitory_conductance, inhibitory_reversal_potential)
            [m, n] = size(obj);
            for i = 1: 1: m
                for j = 1: 1: n
                    obj(i, j).inhibitory_current = inhibitory_conductance.*(obj(i, j).membrane_potential - inhibitory_reversal_potential);
                end
            end
        end

        function obj = compute_polarizations(obj)
            [m, n] = size(obj);
            for i = 1: 1: m
                for j = 1: 1: n
                    delta_v = obj(i, j).membrane_potential - obj(i, j).reference_potential;
                    obj(i, j).depolarization = delta_v.*(delta_v > 0);
                    obj(i, j).hyperpolarization = delta_v.*(delta_v < 0);
                end
            end
        end

    end

end