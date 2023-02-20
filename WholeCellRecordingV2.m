classdef WholeCellRecordingV2
    properties
        %% Acquisition:
        fs

        %% Database:
        filename
        paradigm
        response_duration
        response_samples

        %% Stimulus:
        rate
        npulses

        %% Recording:
        times
        stimulus
        membrane_potential
        injected_current

        %% Cell:
        membrane_capacitance
        input_resistance
        resting_potential
        excitatory_reversal_potential 
        inhibitory_reversal_potential
        threshold_potential
        activation_potential
        steady_state_potential
        alpha_multiplier
        beta_multiplier
        spikes_per_stimulus
        reference_potential

        %% Analysis:
        membrane_current
        leakage_current
        alpha
        beta
        activation_current
        excitatory_conductance
        inhibitory_conductance
        excitatory_current
        inhibitory_current
        depolarizations
        hyperpolarizations

        %% Filters
        differentiatorfir
        lowpassfir

        %% Stats
    end

    methods(Access=public)

        function app = WholeCellRecordingV2(filename, paradigms, response_durations, filter_parameters)
            if nargin > 0
                [m, n] = size(paradigms);
                app(m, n) = app;
                for i = 1: m
                    for j = 1: n
                        app(i, j).filename = filename;
                        app(i, j).paradigm = paradigms(i, j);
                        app(i, j).response_duration = response_durations(i, j);
                    end
                end
                app = app.read_data();
                app = app.read_parameters();
                app = app.build_differentiatorfir(filter_parameters.differentiatorfir);
                app = app.build_lowpassfir(filter_parameters.lowpassfir);
            end
        end

        function app = read_data(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    data = readtable(app(i, j).filename, 'Sheet', app(i, j).paradigm, 'ReadVariableNames', 1, 'VariableNamingRule', 'preserve');
                    try
                        app(i, j).times = data.times;
                        app(i, j).fs = 1/(app(i, j).times(2, 1) - app(i, j).times(1, 1));
                        data.times = [];
                    catch
                        error(strcat('times array was not found in ', app(i, j).paradigm, '!'));
                    end
                    try
                        app(i, j).stimulus = data.stimulus;
                        data.stimulus = [];
                    catch
                        warning(strcat('stimulus array was not found in ', app(i, j).paradigm, '!'));
                        app(i, j).stimulus = [];
                    end
                    try
                        app(i, j).membrane_potential = table2array(data).*1e-3;
                    catch
                        error(strcat('membrane potentials was not found in ', app(i, j).paradigm, '!'));
                    end
                    app(i, j).response_samples = app(i, j).fs * app(i, j).response_duration;
                    if app(i, j).response_samples > size(app(i, j).membrane_potential, 1)
                        app(i, j).response_samples = size(app(i, j).membrane_potential, 1);
                    end
                end
            end
        end

        function app = read_parameters(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    parameters = readtable(app(i, j).filename, 'Sheet', strcat('parameters_', app(i, j).paradigm));
                    app(i, j).injected_current = parameters.Iinj' + zeros(size(app(i, j).membrane_potential));
                    app(i, j).membrane_capacitance = parameters.Cm' + zeros(size(app(i, j).membrane_potential));
                    app(i, j).input_resistance = parameters.Rin' + zeros(size(app(i, j).membrane_potential));
                    app(i, j).resting_potential = parameters.Er' + zeros(size(app(i, j).membrane_potential));
                    app(i, j).excitatory_reversal_potential = parameters.Ee' + zeros(size(app(i, j).membrane_potential));
                    app(i, j).inhibitory_reversal_potential = parameters.Ei' + zeros(size(app(i, j).membrane_potential));
                    app(i, j).threshold_potential = parameters.Et' + zeros(size(app(i, j).membrane_potential));
                    app(i, j).activation_potential = parameters.Eact' + zeros(size(app(i, j).membrane_potential));
                    app(i, j).steady_state_potential = parameters.Ess' + zeros(size(app(i, j).membrane_potential));
                    app(i, j).alpha_multiplier = parameters.xalpha' + zeros(size(app(i, j).membrane_potential));
                    app(i, j).beta_multiplier = parameters.xbeta' + zeros(size(app(i, j).membrane_potential));
                    app(i, j).spikes_per_stimulus = parameters.sps' + zeros(size(app(i, j).membrane_potential));
                    app(i, j).reference_potential = parameters.Eref' + zeros(size(app(i, j).membrane_potential));
                    app(i, j).rate = parameters.rate';
                    app(i, j).npulses = parameters.npulses';
                end
            end
        end

        function app = build_analysis_arrays(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    app(i, j).membrane_current = zeros(size(app(i, j).membrane_potential));
                    app(i, j).leakage_current = zeros(size(app(i, j).membrane_potential));
                    app(i, j).alpha = zeros(size(app(i, j).membrane_potential));
                    app(i, j).beta = zeros(size(app(i, j).membrane_potential));
                    app(i, j).activation_current = zeros(size(app(i, j).membrane_potential));
                    app(i, j).excitatory_conductance = zeros(size(app(i, j).membrane_potential));
                    app(i, j).inhibitory_conductance = zeros(size(app(i, j).membrane_potential));
                    app(i, j).excitatory_current = zeros(size(app(i, j).membrane_potential));
                    app(i, j).inhibitory_current = zeros(size(app(i, j).membrane_potential));
                    app(i, j).depolarizations = zeros(size(app(i, j).membrane_potential));
                    app(i, j).hyperpolarizations = zeros(size(app(i, j).membrane_potential));
                end
            end
        end

        function app = build_differentiatorfir(app, parameters)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    app(i, j).differentiatorfir.window = designfilt('differentiatorfir', ...
                       'FilterOrder', parameters.order, ...
                       'PassbandFrequency', parameters.Fpass, ...
                       'StopbandFrequency', parameters.Fstop, ...
                       'SampleRate', app(i, j).fs);
                    app(i, j).differentiatorfir.delay = mean(grpdelay(app(i, j).differentiatorfir.window));
                end
            end
        end

        function app = build_lowpassfir(app, parameters)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    app(i, j).lowpassfir.window = designfilt('lowpassfir', ...
                        'FilterOrder', parameters.order, ...
                        'PassbandFrequency', parameters.Fpass, ...
                        'StopbandFrequency', parameters.Fstop, ...
                        'SampleRate', app(i, j).fs);
                    app(i, j).lowpassfir.delay = mean(grpdelay(app(i, j).lowpassfir.window));
                end
            end
        end

        function app = update_lowpassfir_filter(app, filter_parameters)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    app(i, j)
                end
            end
        end
    end

    methods
        function app = filter_membrane_potential(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    appender = zeros(app(i, j).lowpassfir.delay*2, size(app(i, j).membrane_potential, 2));
                    Vm = cat(1, appender+app(i, j).membrane_potential(1, :), app(i, j).membrane_potential, appender+app(i, j).membrane_potential(end, :));
                    Vm = filter(app(i, j).lowpassfir.window, Vm);
                    app(i, j).membrane_potential = Vm(app(i, j).lowpassfir.delay*3 + 1: end - app(i, j).lowpassfir.delay, :);
                end
            end
        end

        function app = compute_membrane_current(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    appender = zeros(app(i, j).differentiatorfir.delay*2, ...
                        size(app(i, j).membrane_potential, 2));
                    Vm = cat(1, appender+app(i, j).membrane_potential(1, :), app(i, j).membrane_potential, appender+app(i, j).membrane_potential(end, :));
                    dVmdt = filter(app(i, j).differentiatorfir.window, Vm);
                    app(i, j).membrane_current = app(i, j).membrane_capacitance .* dVmdt(app(i, j).differentiatorfir.delay*3 + 1: end-app(i, j).differentiatorfir.delay, :);
                end
            end
        end

        function app = filter_membrane_current(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    appender = zeros(app(i, j).lowpassfir.delay*2, size(app(i, j).membrane_current, 2));
                    Im = cat(1, appender+app(i, j).membrane_current(1, :), app(i, j).membrane_current, appender+app(i, j).membrane_current(end, :));
                    Im = filter(app(i, j).lowpassfir.window, Im);
                    app(i, j).membrane_current = Im(app(i, j).lowpassfir.delay*3 + 1: end - app(i, j).lowpassfir.delay, :);
                end
            end
        end

    end

    methods
        function app = call(app)
            app = app.filter_membrane_potential();
            app = app.compute_membrane_current();
            app = app.filter_membrane_current();
        end
    end
end