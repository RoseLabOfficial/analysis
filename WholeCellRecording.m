classdef WholeCellRecording
    properties
        %% Acquisition:
        fs

        %% Database:
        filename
        paradigms
        response_durations
        response_samples

        %% Stimulus:
        rates
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
        excitation
        inhibition
        resultant_excitaiton
        resultant_inhibition

        %% Filters
        lowpass
        notch

        %% Stats
        means
        maxima
        minima

        %% Plots
        estimations_figure
        stats_figure
    end

    methods(Access=public)
        %% Datastructure building methods
        function app = WholeCellRecording(filename, paradigms, response_durations, filter_parameters)
            if nargin > 0
                [m, n] = size(paradigms);
                app(m, n) = app;
                for i = 1: m
                    for j = 1: n
                        app(i, j).filename = filename;
                        app(i, j).paradigms = paradigms(i, j);
                        app(i, j).response_durations = response_durations(i, j);
                    end
                end
                app = app.read_data();
                app = app.read_parameters();
                app = app.build_filters(filter_parameters);
                app = app.build_analysis_arrays();
            end
        end

        function app = read_data(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    data = readtable(app(i, j).filename, 'Sheet', app(i, j).paradigms, 'ReadVariableNames', 1, 'VariableNamingRule', 'preserve');
                    try
                        app(i, j).times = data.times;
                        app(i, j).fs = 1/(app(i, j).times(2, 1) - app(i, j).times(1, 1));
                        data.times = [];
                    catch
                        error(strcat('times array was not found in ', app(i, j).paradigms, '!'));
                    end
                    try
                        app(i, j).stimulus = data.stimulus;
                        data.stimulus = [];
                    catch
                        warning(strcat("stimulus array was not found in ", app(i, j).paradigms, "!"));
                        app(i, j).stimulus = [];
                    end
                    try
                        app(i, j).membrane_potential = table2array(data).*1e-3;
                    catch
                        error(strcat('membrane potentials was not found in ', app(i, j).paradigms, '!'));
                    end
                    app(i, j).response_samples = app(i, j).fs * app(i, j).response_durations;
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
                    parameters = readtable(app(i, j).filename, 'Sheet', strcat('parameters_', app(i, j).paradigms));
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
                    app(i, j).rates = parameters.rate';
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
                    app(i, j).excitatory_conductance = zeros(size(app(i, j).membrane_potential, 1), 1);
                    app(i, j).inhibitory_conductance = zeros(size(app(i, j).membrane_potential, 1), 1);
                    app(i, j).excitatory_current = zeros(size(app(i, j).membrane_potential));
                    app(i, j).inhibitory_current = zeros(size(app(i, j).membrane_potential));
                    app(i, j).depolarizations = zeros(size(app(i, j).membrane_potential));
                    app(i, j).hyperpolarizations = zeros(size(app(i, j).membrane_potential));
                    app(i, j).excitation = zeros(size(app(i, j).membrane_potential, 1), 1);
                    app(i, j).inhibition = zeros(size(app(i, j).membrane_potential, 1), 1);
                    app(i, j).resultant_excitaiton = zeros(size(app(i, j).membrane_potential, 1), 1);
                    app(i, j).resultant_inhibition = zeros(size(app(i, j).membrane_potential, 1), 1);
                end
            end
        end

        function app = adjust_membrane_potential_wrt_steady_state(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    app(i, j).injected_current = (1./app(i, j).input_resistance).*(app(i, j).steady_state_potential - app(i, j).resting_potential);
                    app(i, j).membrane_potential = app(i, j).membrane_potential - app(i, j).membrane_potential(1, :) + app(i, j).steady_state_potential;
                end
            end
        end

    end
    %% Filter methods
    methods
        function app = build_notch(app, params)
            nyquist_frequency = app(1, 1).fs/2;
            normalized_stop_frequency = params.Fst/nyquist_frequency;
            normalized_stopband_width = params.BW/nyquist_frequency;
            [app(1, 1).notch.num, app(1, 1).notch.den] = iirnotch(normalized_stop_frequency, normalized_stopband_width);
        end
        
        function app = build_lowpass(app, params)
            app(1, 1).lowpass = designfilt("lowpassiir", ...
                "PassbandFrequency", params.Fp, ...
                "StopbandFrequency", params.Fst, ...
                "PassbandRipple", params.Ap, ...
                "StopbandAttenuation", params.Ast, ...
                "SampleRate", app(1, 1).fs);
        end

        function app = build_filters(app, params)
            app = app.build_lowpass(params.lowpass);
            app = app.build_notch(params.notch);
        end

        function app = filter_membrane_potential(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    appender = zeros(size(app(i, j).membrane_potential));
                    vm = cat(1, appender+app(i, j).membrane_potential(1, :), app(i, j).membrane_potential, appender+app(i, j).membrane_potential(end, :));
                    vm = filtfilt(app(1, 1).lowpass, vm);
                    vm = filtfilt(app(1, 1).notch.num, app(1, 1).notch.den, vm);
                    app(i, j).membrane_potential = vm(size(app(i, j).membrane_potential, 1)+1: end-size(app(i, j).membrane_potential, 1), :);
                end
            end
        end
    end

    %% Computing current methods
    methods

        function app = compute_membrane_current(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    appender = zeros(1, size(app(i, j).membrane_potential, 2));
                      dvmdt = diff(cat(1, appender+app(i, j).membrane_potential(1, :), app(i, j).membrane_potential, appender+app(i, j).membrane_potential(end, :)), 1);
                      dvmdt = dvmdt(2:end, :); 
                      app(i, j).membrane_current = app(i, j).membrane_capacitance.*dvmdt.*app(i, j).fs;
                      app(i, j).membrane_current = app(i, j).membrane_current - app(i, j).membrane_current(1, :);
                end
            end
        end

        function app = compute_leakage_current(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    app(i, j).leakage_current = (1./app(i, j).input_resistance).*(app(i, j).membrane_potential - app(i, j).resting_potential);
                end
            end
        end

        function app = compute_active_conductances(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    app(i, j).alpha = ((1./app(i, j).input_resistance)./(2.*(app(i, j).activation_potential-app(i, j).steady_state_potential)));
                    app(i, j).beta = app(i, j).alpha.*(app(i, j).threshold_potential - app(i, j).steady_state_potential);
                    app(i, j).alpha = app(i, j).alpha.*app(i, j).alpha_multiplier;
                    app(i, j).beta = app(i, j).beta.*app(i, j).beta_multiplier;
                end
            end
        end

        function app = compute_activation_currents(app)
            app = app.compute_active_conductances();
            [m, n] = size(app);
            for i = 1: 1: m
               for j = 1: 1: n
                  app(i, j).activation_current = (app(i, j).alpha.*(app(i, j).membrane_potential-app(i, j).threshold_potential).*(app(i, j).membrane_potential-app(i, j).resting_potential)) ...
                          + (app(i, j).beta.*(app(i, j).membrane_potential-app(i, j).resting_potential));
               end
            end
        end

        function app = compute_passive_conductances(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    samples = size(app(i, j).membrane_potential, 1);
                    A = zeros(2, 2, samples);
                    B = zeros(2, 1, samples);
                    A(1, 1, :) = sum((app(i, j).membrane_potential - app(i, j).excitatory_reversal_potential).^2, 2);
                    A(1, 2, :) = sum((app(i, j).membrane_potential - app(i, j).excitatory_reversal_potential).*(app(i, j).membrane_potential - app(i, j).inhibitory_reversal_potential), 2);
                    A(2, 1, :) = A(1, 2, :);
                    A(2, 2, :) = sum((app(i, j).membrane_potential - app(i, j).inhibitory_reversal_potential).^2, 2);
                    C = app(i, j).membrane_current - app(i, j).injected_current - app(i, j).activation_current + app(i, j).leakage_current;
                    B(1, 1, :) = -sum(C.*(app(i, j).membrane_potential - app(i, j).excitatory_reversal_potential), 2);
                    B(2, 1, :) = -sum(C.*(app(i, j).membrane_potential - app(i, j).inhibitory_reversal_potential), 2);
                    G = pagemtimes(pageinv(A), B);
                    app(i, j).excitatory_conductance = reshape(G(1, 1, :), [samples, 1]);
                    app(i, j).inhibitory_conductance = reshape(G(2, 1, :), [samples, 1]);
                end
            end
        end

        function app = compute_passive_currents(app)
            app = app.compute_passive_conductances();
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    app(i, j).excitatory_current = -1.*app(i, j).excitatory_conductance.*(app(i, j).membrane_potential - app(i, j).excitatory_conductance);
                    app(i, j).inhibitory_current = -1.*app(i, j).inhibitory_conductance.*(app(i, j).membrane_potential - app(i, j).inhibitory_conductance);
                end
            end
        end
    end

    methods

        function app = compute_resultants(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    app(i, j).excitation = (app(i, j).excitatory_conductance > 0).*app(i, j).excitatory_conductance;
                    app(i, j).inhibition = (app(i, j).inhibitory_conductance > 0).*app(i, j).inhibitory_conductance;
                    resultant = app(i, j).excitation - app(i, j).inhibition;
                    app(i, j).resultant_excitaiton = (resultant > 0).*resultant;
                    app(i, j).resultant_inhibition = (resultant < 0).*resultant;
                end
            end
        end

        function app = compute_polarizations(app)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    polarization = app(i, j).membrane_potential - app(i, j).steady_state_potential;
                    app(i, j).depolarizations = (polarization > 0).*polarization;
                    app(i, j).hyperpolarizations = abs((polarization < 0).*polarization);
                end
            end
        end

        function app = compute_means(app, idx)
            [m, n] = size(app);
            for i = 1: m
                for j = 1: n
                    app(i, j).means.depolarization = mean(app(i, j).depolarizations(:, idx), 1);
                    app(i, j).means.hyperpolarization = mean(app(i, j).hyperpolarizations(:, idx), 1);
                    app(i, j).means.excitatory_conductance = mean(app(i, j).excitatory_conductance, 1);
                    app(i, j).means.inhibitory_conductance = mean(app(i, j).inhibitory_conductance, 1);
                    app(i, j).means.excitatory_current = mean(app(i, j).excitatory_current(:, idx), 1);
                    app(i, j).means.inhibitory_current = mean(app(i, j).inhibitory_current(:, idx), 1);
                    app(i, j).means.excitation = mean(app(i, j).excitation, 1);
                    app(i, j).means.inhibition = mean(app(i, j).inhibition, 1);
                    app(i, j).means.resultant_excitation = mean(app(i, j).resultant_excitaiton, 1);
                    app(i, j).means.resultant_inhibition = mean(app(i, j).resultant_inhibition, 1);
                    app(i, j).means.activation_current = mean(app(i, j).activation_current(:, idx), 1);
                    app(i, j).means.membrane_current = mean(app(i, j).membrane_current(:, idx), 1);
                    app(i, j).means.rates = app(i, j).rates(idx, 1);
                    app(i, j).means.npulses = app(i, j).npulses(idx, 1);
                    app(i, j).means.spikes_per_stimulus = app(i, j).spikes_per_stimulus(idx, 1);
                end
            end
        end

        function mean_stats = compute_mean_stats(app)
            means = [app.means];
            pps = [means.rates];
            pulses = [means.npulses];
            sps = [means.spikes_per_stimulus];
            depol = [means.depolarization];
            hyperpol = [means.hyperpolarization];
            ge = [means.excitation];
            gi = [means.inhibition];
            net_ge = [means.resultant_excitation];
            net_gi = [means.resultant_inhibition];
            Ie = [means.excitatory_current];
            Ii = [means.inhibitory_current];
            Iact = [means.activation_current];
            Im = [means.membrane_current];
            mean_stats = [pps; pulses; sps; depol; hyperpol; ge; gi; net_ge; net_gi; Ie; Ii; Iact; Im];
        end

        function meta_stats = compute_meta_stats(app)
            paradigms = [app.paradigms];
            row_names = ["rate"; "pulses"; "sps"; "depol"; "hyperpol"; "ge"; "gi"; "net_ge"; "net_gi"; "Ie"; "Ii"; "Iact"; "Im"];
            meta_stats = array2table(app.compute_mean_stats(), "RowNames", row_names, "VariableNames",paradigms);
        end

    end
    
    %% Call method for estimator
    methods
        function app = call(app)
            app = app.filter_membrane_potential();
            app = app.adjust_membrane_potential_wrt_steady_state();
            app = app.compute_membrane_current();
            app = app.compute_leakage_current();
            app = app.compute_activation_currents();
            app = app.compute_passive_currents();
            app = app.compute_resultants();
            app = app.compute_polarizations();
        end
    end

    %% Plot methods
    methods
        function app = build_plots(app)
            app(1, 1).estimations_figure = figure('Name', strcat('Estimations for ', app(1, 1).filename));
            app(1, 1).stats_figure = figure('Name', strcat('Means for ', app(1, 1).filename));
        end

        function app = plot_estimations(app)
            nplots = 8;
            figure(app(1, 1).estimations_figure);
            [m, n] = size(app);
            for i = 1: m
                tiledlayout(nplots, n);
                ax = cell(nplots, n);
                for k = 1: nplots
                    for j = 1: n
                        ax{j, k} = nexttile;
                        switch k
                            case 1
                                plot(app(i, j).times, app(i, j).membrane_potential);
                                ax{j, k}.Title.String = app(i, j).paradigms;
                                ax{j, k}.Subtitle.String = [strcat("Eact=", num2str(app(i, j).activation_potential(1, :))), strcat('Diff: O=', num2str(size(app(1, 1).lowpass.Coefficients, 1)-1), '; Band=(', num2str(app(1, 1).lowpass.PassbandFrequency), ', ', num2str(app(1, 1).lowpass.StopbandFrequency), ')')];
                                if j == 1
                                    ax{j, k}.YLabel.String = 'Vm (V)';
                                end
                            case 2
                                plot(app(i, j).times, app(i, j).activation_current);
                                ax{j, k}.Subtitle.String = strcat('xalpha=', num2str(app(i, j).alpha_multiplier(1, :)), '; xbeta=', num2str(app(i, j).beta_multiplier(1, :)));
                                if j == 1
                                    ax{j, k}.YLabel.String = 'Iact (A)';
                                end
                            case 3
                                plot(app(i, j).times, app(i, j).leakage_current);
                                ax{j, k}.Subtitle.String = strcat('Rin=', num2str(app(i, j).input_resistance(1, 1)));
                                if j == 1
                                    ax{j, k}.YLabel.String = 'Ileak (A)';
                                end
                            case 4
                                plot(app(i, j).times, app(i, j).membrane_current);
                                ax{j, k}.Subtitle.String = strcat('Cm=', num2str(app(i, j).membrane_capacitance(1, 1)));
                                if j == 1
                                    ax{j, k}.YLabel.String = 'Im (A)';
                                end
                            case 5
                                plot(app(i, j).times, app(i, j).excitatory_conductance, 'r', app(i, j).times, app(i, j).inhibitory_conductance, 'b');
                                ax{j, k}.Subtitle.String = strcat('Ee=', num2str(app(i , j).excitatory_reversal_potential(1, 1)), '; Ei=', num2str(app(i, j).inhibitory_reversal_potential(1, 1)));
                                if j == 1
                                    ax{j, k}.YLabel.String = 'G (S)';
                                end
                            case 6
                                area(app(i, j).times, app(i, j).resultant_excitaiton, 'FaceColor', 'r');
                                hold on;
                                area(app(i, j).times, app(i, j).resultant_inhibition, 'FaceColor', 'b');
                                hold off;
                                if j == 1
                                   ax{j, k}.YLabel.String = 'Resultant G(S)';
                                end
                            case 7
                                plot(app(i, j).times, app(i, j).excitatory_current(:, 1), 'r', app(i, j).times, app(i, j).inhibitory_current(:, 1), 'b');
                                if j == 1
                                    ax{j, k}.YLabel.String = 'Isyn (A)';
                                end
                            case 8
                                if ~(isempty(app(i, j).stimulus))
                                    plot(app(i, j).times, app(i, j).stimulus, 'k');
                                end
                                if j == 1
                                   ax{j, k}.YLabel.String = 'Stimulus (V)';
                                end
                        end
                        linkaxes([ax{j, :}], 'x');
                    end
                    linkaxes([ax{:, k}], 'y');
                end
            end
        end

        function app = plot_meta_stats(app, meta_stats)
            x = categorical(meta_stats.Properties.VariableNames);
            meta_stats = rows2vars(meta_stats);
            nplots = 5;
            figure(app(1, 1).stats_figure);
            tiledlayout(nplots, 1);
            ax = cell(nplots, 1);
            for k = 1: nplots
                ax{k, 1} = nexttile;
                switch k
                    case 1
                        bar(meta_stats.depol, "BarWidth", 0.5, "FaceColor", [1, 0, 0]);
                        hold on;
                        bar(-1.*meta_stats.hyperpol, "BarWidth", 0.5, "FaceColor", [0, 0, 1]);
                        ylim([-1*max(max(meta_stats.depol), max(meta_stats.hyperpol)),max(max(meta_stats.depol), max(meta_stats.hyperpol))]);
                        hold off;
                        ylabel("polarization (V)")
                    case 2
                        bar(meta_stats.ge, "BarWidth", 0.5, "FaceColor", [1, 0, 0]);
                        hold on;
                        bar(-1.*meta_stats.gi, "BarWidth", 0.5, "FaceColor", [0, 0, 1]);
                        bar(meta_stats.net_ge, "BarWidth", 0.25, "FaceColor", [1, 0.5, 0]);
                        bar(meta_stats.net_gi, "BarWidth", 0.25, "FaceColor", [0, 0.5, 1]);
                        ylim([-1*max(max(meta_stats.ge), max(meta_stats.gi)),max(max(meta_stats.ge), max(meta_stats.gi))]);
                        hold off;
                        ylabel("syn conductances (S)")
                    case 3
                        bar(meta_stats.Ie, "BarWidth", 0.5, "FaceColor", [1, 0, 0]);
                        hold on;
                        bar(-1.*meta_stats.Ii, "BarWidth", 0.5, "FaceColor", [0, 0, 1]);
                        ylim([-1*max(max(meta_stats.Ie), max(meta_stats.Ii)),max(max(meta_stats.Ie), max(meta_stats.Ii))]);
                        hold off;
                        ylabel("syn currents (A)")
                    case 4
                        bar(meta_stats.Iact, "BarWidth", 0.5, "FaceColor", [1, 0, 0]);                      
                        ylim([-1*max(meta_stats.Iact),max(meta_stats.Iact)]);
                        ylabel("act current (A)")
                    case 5
                        bar(meta_stats.Im, "BarWidth", 0.5, "FaceColor", [1, 0, 0]);                      
                        ylim([-1*max(meta_stats.Im),max(meta_stats.Im)]);
                        ylabel("mem currents (A)")
                end
                xticklabels(x);
                yyaxis right;
                plot(meta_stats.sps, '-ok', 'MarkerFaceColor',[0, 0, 0], 'MarkerSize', 6);
                ylim([-1*(max(meta_stats.sps)+0.2), (max(meta_stats.sps)+0.2)]);
                ylabel("sps");
            end
        end

    end
end