classdef WholeCellRecording
   properties
       %% Sampling/Acquistion properties.
       Fs % Sampling Frequency of acquisition (Hz).
       response_samples % How long the response lasts in a sampling domain (samples).

       %% Stimulus properties.
       rate
       npulses
       
       %% Cell properties.
       Cm % Membrane Capacitance (F).
       Rin % Input Resistance of the cell (Ohms).
       Er % Resting Potential of the cell: Steady State Potential at 0nA current clamp(V).
       Ee % Excitatory Reversal Potential: potential at which excitation is reversed and appears as hyperpolarization (V).
       Ei % Inhibitory Reversal Potential: potential at which inhibition is reversed and appears as depolarization (V).
       Et % Threshold Potential: potential at which an action potential is triggered or supra-threshold active currents are triggered (V).
       Eact % Activation Potential: potential at which subthreshold active currents are triggered (V).
       Ess % Steady State Potential: settled potential at a given current clamp when no stimulus is presented (V).
       xalpha % Alpha Multiplier: multiplier on the active current estimate (dimensionless).
       xbeta % Beta Multiplier: multiplier on the compensator for the active current (dimensionless).
       sps % spikes for stimulus repetition (spikes/stim.rep.).
       Eref % Reference Potential: potential used as reference to measure de- and hyperpolarizations (V).
       
       %% filter properties
       CutOffFrequency % Frequency in the signal where the noise starts to dominate the signal (Hz). 
       FilterOrder % The number of previous timesteps that need to be used in filtering (dimensionless).
       PassbandRipple % The amplitude of ripples in frequencies of the passband (dB).
       StopbandAttenuation % The attentuation of noise (dB).
       
       %% analysis properties
       stimulus
       time % times(secs) extracted from recorded data.
       Vm % Transmembrane Potential (V) recorded at Fs.
       Im % Transmembrane current (A) computed as Cm*(dVm/dt).
       Iinj % Injected Current (A).
       Ileak % Leakage current (A) of the cell, computed as (1/Rin)*(Vm - Er).
       alpha % Active conductance constant (S/V); Estimated using algorithm from Alluri, 2021. 
       beta % Active conductance constant (S); Estimated using algorithm from Alluri, 2021.
       Iactive % Active current (A) computed using active conductance terms. Alluri, 2021.
       ge % Excitatory conductance (S) estimated using alogrithm from Alluri, 2016.
       gi % Inhibitory conductance (S) estimated using alogrith from Alluri, 2016.
       Ie % Excitatory currents (A) at various current clamps computed as ge*(Vm - Ee).
       Ii % Inhibitory currents (A) at various current clamps computed as gi*(Vm - Ei).
       depolarizations % Vm > Er at every time point (V).
       hyperpolarizations % Vm < Er at every time point (V).
       excitation % ge > 0. Biological data is noisy, excitation variable is free of aberrent negative values of ge.
       inhibition % gi > 0. Biological data is noisy, inhibition variable is free of aberrent negative values of gi after non-linear filtering using alpha and beta.
       paradigm % paradigm (string) is different types of stimuli, for pulse rate 5pps, 10pps, 60pps, etc., or for duration 20ms, 40ms, 160ms, etc.
       filename % filename (string) is the name of the file containing the current clamp data along with various parameters for cell and membrane potential constants.
       response_durations

       %% stats
        ge_net
        gi_net
        ge_mean
        gi_mean
        depolarizations_mean
        hyperpolarizations_mean
        Iactive_mean
        sps_mean

        %% meta stats
        maxima
        minima
        drop_3dB
   end
   properties(Access=private)
       %% fig properties
       fig
   end
   %% Methods to access data.
   methods(Access=private)
       function app = readXLSXdata(app)
           [m, n] = size(app);
           for i = 1: 1: m
               for j = 1: 1: n
                    data = readtable(app(i, j).filename, 'Sheet', app(i, j).paradigm, 'ReadVariableNames', 1, 'VariableNamingRule', 'preserve');
                    try
                        app(i, j).time = data.times;
                        data.times = [];
                    catch
                        error(strcat('times array was not found! Check ', app(i, j).filename, ' and ', app(i, j).paradigm, ' !'));
                    end
                    try
                        app(i, j).stimulus = data.stimulus;
                        data.stimulus = [];
                    catch
                        app(i, j).stimulus = zeros(size(app(i, j).time));
                    end
                    app(i, j).Vm = table2array(data).*1e-3;
                    app(i, j).Fs = 1/(app(i, j).time(2, 1) - app(i, j).time(1, 1));
                    app(i, j).response_samples = floor(app(i, j).Fs*app(i, j).response_durations(i));
                    samples = size(app(i, j).time, 1);
                    if app(i, j).response_samples > samples-1
                       app(i, j).response_samples = samples-1;
                    end
               end
           end
       end

       function app = adjust_membrane_potential_with_steady_state(app)
           [m, n] = size(app);
           for i = 1: 1: m
               for j = 1: 1: n
                   app(i, j).Iinj = (1./app(i, j).Rin).*(app(i, j).Ess - app(i, j).Er);
                   app(i, j).Vm = app(i, j).Vm - app(i, j).Vm(1, :) + app(i, j).Ess;
               end
           end
       end

       function app = readXLSXparameters(app)
           [m, n] = size(app);
           for i = 1: 1: m
               for j = 1: 1: n
                    parameters = readtable(app(i, j).filename, 'Sheet', strcat("parameters_", app(i, j).paradigm));
                    app(i, j).Iinj = parameters.Iinj' + zeros(size(app(i, j).Vm));
                    app(i, j).Cm = parameters.Cm' + zeros(size(app(i, j).Vm));
                    app(i, j).Rin = parameters.Rin' + zeros(size(app(i, j).Vm));
                    app(i, j).Er = parameters.Er' + zeros(size(app(i, j).Vm));
                    app(i, j).Ee = parameters.Ee' + zeros(size(app(i, j).Vm));
                    app(i, j).Ei = parameters.Ei' + zeros(size(app(i, j).Vm));
                    app(i, j).Et = parameters.Et' + zeros(size(app(i, j).Vm));
                    app(i, j).Eact = parameters.Eact' + zeros(size(app(i, j).Vm));
                    app(i, j).Ess = parameters.Ess' + zeros(size(app(i, j).Vm));
                    app(i, j).xalpha = parameters.xalpha' + zeros(size(app(i, j).Vm));
                    app(i, j).xbeta = parameters.xbeta' + zeros(size(app(i, j).Vm));
                    app(i, j).sps = parameters.sps' + zeros(size(app(i, j).Vm));
                    app(i, j).Eref = parameters.Eref' + zeros(size(app(i, j).Vm));
                    app(i, j).rate = parameters.rate(1);
                    app(i, j).npulses = parameters.npulses(1);
               end
           end
       end

       function app = build_arrays(app)
           [m, n] = size(app);
           for i = 1: 1: m
               for j = 1: 1: n
                    samples = size(app(i, j).Vm, 2);
                    app(i, j).Iactive = zeros(size(app(i, j).Vm));
                    app(i, j).ge = zeros(samples, 1);
                    app(i, j).gi = zeros(samples, 1);
                    app(i, j).Ie = zeros(size(app(i, j).Vm));
                    app(i, j).Ii = zeros(size(app(i, j).Vm));
                    app(i, j).depolarizations = zeros(samples, 1);
                    app(i, j).hyperpolarizations = zeros(samples, 1);
                    app(i, j).excitation = zeros(samples, 1);
                    app(i, j).inhibition = zeros(samples, 1);
                    app(i, j).ge_mean = 0;
                    app(i, j).gi_mean = 0;
                    app(i, j).ge_net = 0;
                    app(i, j).gi_net = 0;
                    app(i, j).Iactive_mean = 0;
                    app(i, j).depolarizations_mean = 0;
                    app(i, j).hyperpolarizations_mean = 0;
                    app(i, j).sps_mean = 0;
                    app(i, j).maxima.ge_mean = zeros(1, 2);
                    app(i, j).maxima.gi_mean = zeros(1, 2);
                    app(i, j).maxima.ge_net = zeros(1, 2);
                    app(i, j).maxima.gi_net = zeros(1, 2);
                    app(i, j).maxima.Iactive_mean = zeros(1, 2);
                    app(i, j).maxima.depolarizations_mean = zeros(1, 2);
                    app(i, j).maxima.hyperpolarizations_mean = zeros(1, 2);
                    app(i, j).maxima.sps_mean = zeros(1, 2);
                    app(i, j).minima.ge_mean = zeros(1, 2);
                    app(i, j).minima.gi_mean = zeros(1, 2);
                    app(i, j).minima.ge_net = zeros(1, 2);
                    app(i, j).minima.gi_net = zeros(1, 2);
                    app(i, j).minima.Iactive_mean = zeros(1, 2);
                    app(i, j).minima.depolarizations_mean = zeros(1, 2);
                    app(i, j).minima.hyperpolarizations_mean = zeros(1, 2);
                    app(i, j).minima.sps_mean = zeros(1, 2);
                    app(i, j).drop_3dB.ge_mean = zeros(1, 2);
                    app(i, j).drop_3dB.gi_mean = zeros(1, 2);
                    app(i, j).drop_3dB.ge_net = zeros(1, 2);
                    app(i, j).drop_3dB.gi_net = zeros(1, 2);
                    app(i, j).drop_3dB.Iactive_mean = zeros(1, 2);
                    app(i, j).drop_3dB.depolarizations_mean = zeros(1, 2);
                    app(i, j).drop_3dB.hyperpolarizations_mean = zeros(1, 2);
                    app(i, j).drop_3dB.sps_mean = zeros(1, 2);
               end
           end
       end
   end

   %% Methods for reverse estimation:
   methods(Access=public)
       function app = reset_passive_conductances(app, ge_gain, gi_gain)
           tStart = tic;
           [m, n] = size(app);
           for i = 1: 1: m
               for j = 1: 1: n
                   app(i, j).ge = app(i, j).ge - app(i, j).ge(1, :);
                   app(i, j).gi = app(i, j).gi - app(i, j).gi(1, :);
                   app(i, j).ge = ge_gain.*app(i, j).ge;
                   app(i, j).gi = gi_gain.*app(i, j).gi;
               end
           end
           fprintf('[%d secs] Resetting passive conductances using user conductance gains (%d, %d). \n', toc(tStart), ge_gain, gi_gain);
       end

       function app = reverse_estimate_membrane_potential(app, ge_gain, gi_gain)
            app = app.reset_passive_conductances(ge_gain, gi_gain);
            tStart = tic;
            [m, n] = size(app);
            for i = 1: 1: m
                for j = 1: 1: n
                   app(i, j).Vm = 0.*app(i, j).Vm + app(i, j).Ess;
                   for k = 2: 1: size(app(i, j).Vm, 1)
                       app(i, j).Vm(k, :) = app(i, j).Vm(k-1, :) + (1./app(i, j).Fs).*(1./app(i, j).Cm(k-1, :)).*(app(i, j).Iinj(k-1, :) ...
                           - (1./app(i, j).Rin(k-1, :)).*(app(i, j).Vm(k-1, :) - app(i, j).Er(k-1, :)) ...
                           - app(i, j).ge(k-1, :).*(app(i, j).Vm(k-1, :) - app(i, j).Ee(k-1, :)) ...
                           - app(i, j).gi(k-1, :).*(app(i, j).Vm(k-1, :) - app(i, j).Ei(k-1, :)))...
                          - app(i, j).alpha(k-1, :).*(app(i, j).Vm(k-1, :) - app(i, j).Et(k-1, :)).*(app(i, j).Vm(k-1, :) - app(i, j).Er(k-1, :)) ...
                           - app(i, j).beta(k-1, :).*(app(i, j).Vm(k-1, :) - app(i, j).Er(k-1, :));
                   end
                end
            end
            fprintf('[%d secs] Reverse estimating Vm using user conductance gains (%d, %d). \n', toc(tStart), ge_gain, gi_gain);
       end
       function app = plot_reverse_estimations(app)
            tStart = tic;
            nplots = 2;
            [m, n] = size(app);
            app(1, 1).fig = figure('Name', strcat(app(1, 1).filename, ' Reconstructions'));
            for i = 1: 1: m
                tiledlayout(nplots, n);
                ax = cell(nplots, n);
                for k = 1: 1: nplots
                    for j = 1: 1: n
                        ax{j, k} = nexttile;
                        switch k
                            case 1
                                plot(app(i, j).time, app(i, j).Vm);
                                ax{j, k}.Title.String = app(i, j).paradigm;
                                if j == 1
                                    ax{j, k}.YLabel.String = 'Vm (V)';
                                end
                            case 2
                                plot(app(i, j).time, app(i, j).ge, 'r', app(i, j).time, app(i, j).gi, 'b');
                                if j == 1
                                    ax{j, k}.YLabel.String = 'G (S)';
                                end
                                ax{j, k}.XLabel.String = 'time (sec)'; 
                        end
                        linkaxes([ax{j, :}], 'x');
                    end
                end
            end
            fprintf('[%d secs] Plotting reverse estimations. \n', toc(tStart));
        end
   end
   %% Methods to filter data
   methods(Access=public)
       function app = zero_phase_filter_Vm(app, filter_parameters)
           tStart = tic;
           [m, n] = size(app);
           for i = 1: 1: m
               for j = 1: 1: n
                   Fnorm = filter_parameters.CutOffFrequency/(app(i, j).Fs/2);
                   lpFilt = designfilt('lowpassiir', ...
                                'PassbandFrequency', Fnorm, ...
                                'FilterOrder', filter_parameters.FilterOrder, ...
                                'PassbandRipple', filter_parameters.PassbandRipple, ...
                                'StopbandAttenuation', filter_parameters.StopbandAttenuation);
                   v = cat(1, ones(size(app(i, j).Vm)).*(app(i, j).Vm(1, :)), app(i, j).Vm);
                   v = cat(1, v, ones(size(app(i, j).Vm)).*(app(i, j).Vm(end, :)));
                   vn = v - v(1, :);
                   Vmf = filtfilt(lpFilt, vn);
                   app(i, j).Vm = Vmf(size(app(i, j).Vm, 1)+1:end-size(app(i, j).Vm, 1), :) + app(i, j).Vm(1, :);
               end
           end
           fprintf('[%d secs] Zero phase filtering Vm \n', toc(tStart));
       end
       function app = zero_phase_filter_Im(app, filter_parameters)
            tStart = tic;
            [m, n] = size(app);
            for i = 1: 1: m
                for j = 1: 1: n
                    Fnorm = filter_parameters.CutOffFrequency2/(app(i, j).Fs/2);
                    lpFilt = designfilt('lowpassiir', ...
                                'PassbandFrequency', Fnorm, ...
                                'FilterOrder', filter_parameters.FilterOrder, ...
                                'PassbandRipple', filter_parameters.PassbandRipple, ...
                                'StopbandAttenuation', filter_parameters.StopbandAttenuation);
                    I = cat(1, ones(size(app(i, j).Im)).*(app(i, j).Im(1, :)), app(i, j).Im);
                    I = cat(1, I, ones(size(app(i, j).Im)).*(app(i, j).Im(end, :)));
                    in = I - I(1, :);
                    Imf = filtfilt(lpFilt, in);
                    app(i, j).Im = Imf(size(app(i, j).Im, 1)+1:end-size(app(i, j).Im, 1), :) + app(i, j).Im(1, :);
                end
            end
            fprintf('[%d secs] Zero phase filtering Im \n', toc(tStart));
       end
   end
   %% Methods for conductance reconstruction.
   methods(Access=public)
       function app = compute_active_conductances(app)
           tStart = tic;
           [m, n] = size(app);
           for i = 1: 1: m
               for j = 1: 1: n
                   app(i, j).alpha = ((1./app(i, j).Rin)./(2.*(app(i, j).Eact-app(i, j).Ess)));
                   app(i, j).beta = app(i, j).alpha.*(app(i, j).Et - app(i, j).Ess);
                   app(i, j).alpha = app(i, j).alpha.*app(i, j).xalpha;
                   app(i, j).beta = app(i, j).beta.*app(i, j).xbeta;
               end
           end
           fprintf('[%d secs] Computed active conductances (constants)\n', toc(tStart));
       end

       function app = compute_active_currents(app)
           tStart = tic;
           [m, n] = size(app);
           for i = 1: 1: m
               for j = 1: 1: n
                  app(i, j).Iactive = (app(i, j).alpha.*(app(i, j).Vm-app(i, j).Et).*(app(i, j).Vm-app(i, j).Er)) ...
                          + (app(i, j).beta.*(app(i, j).Vm-app(i, j).Er));
               end
           end
           fprintf('[%d secs] Computed active currents\n', toc(tStart));
       end

       function app = compute_leakage_currents(app)
           tStart = tic;
           [m, n] = size(app);
           for i = 1: 1: m
               for j = 1: 1: n
                   app(i, j).Ileak = (1./app(i, j).Rin).*(app(i, j).Vm-app(i, j).Er);
               end
           end
           fprintf('[%d secs] Computed leakage currents\n', toc(tStart));
       end

       function app = compute_membrane_currents(app)
           tStart=tic;
           [m, n] = size(app);
           for i = 1: 1: m
               for j = 1: 1: n
                    Vm_appended = cat(1, app(i, j).Vm(1, :), app(i, j).Vm, app(i, j).Vm(end, :));
                    dVmdt = diff(Vm_appended);
                    app(i, j).Im = app(i, j).Cm.*dVmdt(1:end-1, :).*(app(i, j).Fs);
               end
           end
           fprintf('[%d secs] Computed membrane currents\n', toc(tStart));
       end
       function app = compute_passive_conductances(app)
           tStart = tic;
           [m, n] = size(app);
           for i = 1: 1: m
               for j = 1: 1: n
                   samples = size(app(i, j).Vm, 1);
                   A = zeros(2, 2, samples);
                   B = zeros(2, 1, samples);
                   A(1, 1, :) = sum((app(i, j).Vm - app(i, j).Ee).^2, 2);
                   A(1, 2, :) = sum((app(i, j).Vm - app(i, j).Ee).*(app(i, j).Vm - app(i, j).Ei), 2);
                   A(2, 1, :) = A(1, 2, :);
                   A(2, 2, :) = sum((app(i, j).Vm - app(i, j).Ei).^2, 2);
                   C = app(i, j).Im - app(i, j).Iinj - app(i, j).Iactive + app(i, j).Ileak;
                   B(1, 1, :) = -sum(C.*(app(i, j).Vm - app(i, j).Ee), 2);
                   B(2, 1, :) = -sum(C.*(app(i, j).Vm - app(i, j).Ei), 2);
                   G = pagemtimes(pageinv(A), B);
                   app(i, j).ge = reshape(G(1, 1, :), [samples, 1]);
                   app(i, j).gi = reshape(G(2, 1, :), [samples, 1]);
               end
           end
           fprintf('[%d secs] Computed passive conductances\n', toc(tStart));
       end

       function app = compute_passive_conductances_nopage(app)
           tStart = tic;
           [m, n] = size(app);
           for i = 1: 1: m
               for j = 1: 1: n
                   samples = size(app(i, j).Vm, 1);
                   A = zeros(2, 2, samples);
                   B = zeros(2, 1, samples);
                   G = zeros(2, 1, samples);
                   A(1, 1, :) = sum((app(i, j).Vm - app(i, j).Ee).^2, 2);
                   A(1, 2, :) = sum((app(i, j).Vm - app(i, j).Ee).*(app(i, j).Vm - app(i, j).Ei), 2);
                   A(2, 1, :) = A(1, 2, :);
                   A(2, 2, :) = sum((app(i, j).Vm - app(i, j).Ei).^2, 2);
                   C = app(i, j).Im - app(i, j).Iinj - app(i, j).Iactive + app(i, j).Ileak;
                   B(1, 1, :) = -sum(C.*(app(i, j).Vm - app(i, j).Ee), 2);
                   B(2, 1, :) = -sum(C.*(app(i, j).Vm - app(i, j).Ei), 2);
                   for k = 1: 1: samples
                       G(:, :, k) = pinv(A(:, :, k))*B(:, :, k);
                   end
%                    G = pagemtimes(pageinv(A), B);
                   app(i, j).ge = reshape(G(1, 1, :), [samples, 1]);
                   app(i, j).gi = reshape(G(2, 1, :), [samples, 1]);
               end
           end
           fprintf('[%d secs] Computed passive conductances\n', toc(tStart));
       end

       function app = compute_passive_currents(app)
           tStart = tic;
               [m, n] = size(app);
               for i = 1: 1: m
                   for j = 1: 1: n
                       app(i, j).Ie = app(i, j).ge.*(app(i, j).Vm - app(i, j).Ee);
                       app(i, j).Ii = app(i, j).gi.*(app(i, j).Vm - app(i, j).Ei);
                   end
               end
           fprintf('[%d secs] Computed passive currents\n', toc(tStart));
       end
   end
   %% Methods for stats
   methods(Access=public)
       function app = compute_stats(app)
           tStart = tic;
           [m, n] = size(app);
           for i = 1: 1: m
               for j = 1: 1: n
                   del_Vm = app(i, j).Vm - app(i, j).Eref;
                   app(i, j).depolarizations = del_Vm.*(del_Vm(:, 1) > 0);
                   app(i, j).hyperpolarizations = del_Vm.*(del_Vm(:, 1) < 0);
                   app(i, j).excitation = app(i, j).ge.*(app(i, j).ge>0);
                   app(i, j).inhibition = app(i, j).gi.*(app(i, j).gi>0);
%                    resultant_conductance = app(i, j).excitation(1:app(i, j).response_samples) - app(i, j).inhibition(1:app(i, j).response_samples);
                   resultant_conductance = app(i, j).excitation - app(i, j).inhibition;
                   app(i, j).ge_net = mean(resultant_conductance.*(resultant_conductance>0), 1);
                   app(i, j).gi_net = -1*mean(resultant_conductance.*(resultant_conductance<0), 1);
                   app(i, j).ge_mean = mean(app(i, j).excitation(1:app(i, j).response_samples), 1);
                   app(i, j).gi_mean = mean(app(i, j).inhibition(1:app(i, j).response_samples), 1);
                   app(i, j).Iactive_mean = mean(app(i, j).Iactive(1:app(i, j).response_samples, 1), 1);
                   app(i, j).depolarizations_mean = mean(app(i, j).depolarizations(1:app(i, j).response_samples, 1), 1);
                   app(i, j).hyperpolarizations_mean = mean(app(i, j).hyperpolarizations(1:app(i, j).response_samples, 1), 1);
                   app(i, j).sps_mean = mean(app(i, j).sps(1:app(i, j).response_samples, 1), 1);
               end
           end
           fprintf('[%d secs] Computed Stats\n', toc(tStart));
       end

       function stats = get_stats(app)
           measures = ["rate", "ge_net", "gi_net", "ge_mean", "gi_mean", "depolarizations", "hyperpolarizations", "sps", "Iactive"];
           experiments = [app.paradigm];
           rates = [app.rate];
           ge_nets = [app.ge_net];
           gi_nets = [app.gi_net];
           ge_means = [app.ge_mean];
           gi_means = [app.gi_mean];
           depolarizations_means = [app.depolarizations_mean];
           hyperpolarizations_means = [app.hyperpolarizations_mean];
           sps_means = [app.sps_mean];
           Iactive_means = [app.Iactive_mean];
           stats = [rates; ge_nets; gi_nets; ge_means; gi_means; depolarizations_means; hyperpolarizations_means; sps_means; Iactive_means];
           stats = array2table(stats, "RowNames",measures,"VariableNames",experiments);
       end

       function plot_stats(app)
            rates = categorical(arrayfun(@(a)num2str(a), [app.rate], 'UniformOutput', 0));
            figure();
            hold on;
            bar(rates, [app.ge_mean], "BarWidth", 0.5, "FaceColor", [0.5, 0.0, 0.0]);
            bar(rates, [app.ge_net], "BarWidth", 0.35, "FaceColor", [1.0, 0.0, 0.0]);
            bar(rates, -1.*[app.gi_mean], "BarWidth", 0.5, "FaceColor", [0, 0.0, 0.7]);
            bar(rates, -1.*[app.gi_net], "BarWidth", 0.35, "FaceColor", [0.0, 0.5, 1.0]);
            hold off;
            ylim([-1*max(max([app.ge_mean]), max([app.gi_mean])),max(max([app.ge_mean]), max([app.gi_mean]))]);
            yyaxis right;
            plot(rates, [app.sps_mean], '-ok', 'MarkerFaceColor',[0, 0, 0], 'MarkerSize', 6);
            ylim([-1*(max([app.sps_mean])+0.2), (max([app.sps_mean])+0.2)]);
       end

       function write_stats_to_file(~, stats, filename, SheetName)
            tStart = tic;
            writetable(stats, filename, 'Sheet', SheetName, 'WriteRowNames', 1);
            fprintf('[%d secs] Writing Stats to %s\n', toc(tStart), filename);
        end
   end
   %% Methods for meta stats.
   methods(Access=public)
        function meta_stats = compute_meta_stats(app)
           tStart = tic;
           RowNames = ["ge_net", "gi_net", "ge_mean", "gi_mean", "depolarizations", "hyperpolarizations", "sps", "Iactive"];
           VariableNames = ["min value", "min rate", "max value", "max rate", "Percent change", "-3dB value", "lower cutoff", "upper cutoff", "Q"];
           rates = [app.rate];
           ge_nets = [app.ge_net];
           gi_nets = [app.gi_net];
           ge_means = [app.ge_mean];
           gi_means = [app.gi_mean];
           Iactive_means = [app.Iactive_mean];
           depolarizations_means = [app.depolarizations_mean];
           hyperpolarizations_means = [app.hyperpolarizations_mean];
           sps_means = [app.sps_mean];
           ge_net_points = app.estimate_filter_points(ge_nets, rates);
           gi_net_points = app.estimate_filter_points(gi_nets, rates);
           ge_mean_points = app.estimate_filter_points(ge_means, rates);
           gi_mean_points = app.estimate_filter_points(gi_means, rates);
           depolarizations_mean_points = app.estimate_filter_points(depolarizations_means, rates);
           hyperpolarizations_mean_points = app.estimate_filter_points(-1.*hyperpolarizations_means, rates);
           hyperpolarizations_mean_points(:, 1) = -1.*hyperpolarizations_mean_points(:, 1);
           hyperpolarizations_mean_points(:, 3) = -1.*hyperpolarizations_mean_points(:, 3);
           sps_mean_points = app.estimate_filter_points(sps_means, rates);
           Iactive_mean_points = app.estimate_filter_points(Iactive_means, rates);
           meta_stats = [ge_net_points; gi_net_points; ge_mean_points; gi_mean_points; depolarizations_mean_points; hyperpolarizations_mean_points; sps_mean_points; Iactive_mean_points];
           meta_stats = array2table(meta_stats, "RowNames", RowNames, "VariableNames",VariableNames);
           fprintf('[%d secs] Computed Meta Stats\n', toc(tStart));
        end
        
        function meta_stats = get_meta_stats(app)
            meta_stats = app.compute_meta_stats();
        end

        function write_meta_stats_to_file(~, stats, filename, SheetName)
            tStart = tic;
            writetable(stats, filename, 'Sheet', SheetName, 'WriteRowNames', 1);
            fprintf('[%d secs] Writing Meta Stats to %s\n', toc(tStart), filename);
        end

        function points = estimate_filter_points(~, values, rates)
            [max_value, max_index] = max(values, [], 2);
            [min_value, min_index] = min(values, [], 2);
            percent_change = ((max_value - min_value)/max_value)*100;
            drop_3dB_value = max_value*0.707;
            lower_values = values(1:max_index);
            lower_rates = rates(1:max_index);
            upper_values = values(max_index:end);
            upper_rates = rates(max_index:end);
            if length(lower_values(:)) > 1
                [lower_values, m1, ~] = unique(lower_values, 'first');
                [c1, d1] = sort(m1);
                lower_values = lower_values(d1);
                lower_rates = lower_rates(c1);
                if length(c1(:)) > 1
                    lower_cutoff = interp1(lower_values, lower_rates, drop_3dB_value, "linear", "extrap");
                    if lower_cutoff < 0
                        lower_cutoff = 0;
                    end
                else
                    lower_cutoff = lower_rates;
                end
            else
                lower_cutoff = lower_rates;
            end
            if length(upper_values(:)) > 1
                [upper_values, m2, ~] = unique(upper_values, 'last');
                [c2, d2] = sort(m2);
                upper_values = upper_values(d2);
                upper_rates = upper_rates(c2);
                if length(c2(:)) > 1
                    upper_cutoff = interp1(upper_values, upper_rates, drop_3dB_value, "linear", "extrap");
                else
                    upper_cutoff = upper_rates;
                end
            else
                upper_cutoff = upper_rates;
            end
            Q = rates(max_index)/(2*(upper_cutoff - lower_cutoff));
            points = [min_value, rates(min_index), max_value, rates(max_index), percent_change, drop_3dB_value, lower_cutoff, upper_cutoff, Q];
        end
       
        function cutoff_rate = estimate_cutoff_rate(~, rates, values, query_value)
               high_index = find(values >= query_value, 1, 'first');
               if isempty(high_index)
                   high_index = size(rates, 2);
               end
               low_index = find(values <= query_value, 1, 'last');
               if isempty(low_index)
                   low_index = 1;
               end
               if high_index == low_index
                   cutoff_rate = rates(high_index);
               else
                   cutoff_rate = rates(low_index) + (query_value - values(low_index))*(rates(high_index) - rates(low_index))/(values(high_index) - values(low_index));
               end
           end
        
           function  points = compute_3dB_points(app, values, rates)
               [max_value, max_index] = max(values, [], 2);
               drop_3dB_value = max_value*0.707;
               drop_3dB_rate = app.estimate_cutoff_rate(rates, values, drop_3dB_value);
               Q = rates(max_index)/(2*abs(rates(max_index) - drop_3dB_rate));
               points = [max_value, rates(max_index), drop_3dB_value, drop_3dB_rate, Q];
           end
   end
   %% Methods for plotting.
   methods(Access=public)
        function app = plots(app)
            tStart = tic;
            nplots = 6;
            [m, n] = size(app);
            app(1, 1).fig = figure('Name', strcat(app(1, 1).filename, ' Reconstructions'));
            for i = 1: 1: m
                tiledlayout(nplots, n);
                ax = cell(nplots, n);
                for k = 1: 1: nplots
                    for j = 1: 1: n
                        ax{j, k} = nexttile;
                        switch k
                            case 1
                                plot(app(i, j).time, app(i, j).Vm);
                                ax{j, k}.Title.String = app(i, j).paradigm;
                                if j == 1
                                    ax{j, k}.YLabel.String = 'Vm (V)';
                                end
                            case 2
                                plot(app(i, j).time, app(i, j).Iactive);
                                if j == 1
                                    ax{j, k}.YLabel.String = 'Iactive (A)';
                                end
                            case 3
                                plot(app(i, j).time, app(i, j).Ileak);
                                if j == 1
                                    ax{j, k}.YLabel.String = 'Ileak (A)';
                                end
                            case 4
                                plot(app(i, j).time, app(i, j).Im);
                                if j == 1
                                    ax{j, k}.YLabel.String = 'Im (A)';
                                end
                            case 5
                                plot(app(i, j).time, app(i, j).ge, 'r', app(i, j).time, app(i, j).gi, 'b');
                                if j == 1
                                    ax{j, k}.YLabel.String = 'G (S)';
                                end
                            case 6
                                plot(app(i, j).time, app(i, j).Ie(:, 1), 'r', app(i, j).time, app(i, j).Ii(:, 1), 'b');
                                if j == 1
                                    ax{j, k}.YLabel.String = 'Isyn (A)';
                                end
                                ax{j, k}.XLabel.String = 'time (sec)';
                        end
                        linkaxes([ax{j, :}], 'x');
                    end
                    linkaxes([ax{:, k}], 'y');
                end
            end
            fprintf('[%d secs] Plotting data\n', toc(tStart));
        end
   end
   %% Methods for constructor and run.
   methods(Access=public)
       %% Constructor
       function app = WholeCellRecording(filename, paradigms, response_durations)
            if nargin > 0
                tStart = tic;
                %% Input format check
                if ~isa(filename, 'string')
                    error('filename must be a string'); 
                end
                if ~isa(paradigms, 'string')
                    error('conditions (rates/durations) must be strings'); 
                end
                if ~isa(response_durations, 'double')
                    error('responseDurations must be doubles');
                end
                if ~isequal(size(paradigms), size(response_durations))
                    error('Dimensions of conditions and responseDurations must match.');
                end
                [m, n] = size(response_durations);
                app(m, n) = app;
                for i = 1: 1: m
                    for j = 1: 1: n
                        app(i, j).filename = filename;
                        app(i, j).paradigm = paradigms(i, j);
                        app(i, j).response_durations = response_durations(i, j);
                    end
                end
                %% Reading worksheets from excel files.
                app = app.readXLSXdata();
                app = app.build_arrays();
                app = app.readXLSXparameters();
                app = app.adjust_membrane_potential_with_steady_state();
                fprintf('[%d secs] Read %s\n', toc(tStart), filename);
            end
       end

       function app = call(app, filter_parameters)
            app = app.zero_phase_filter_Vm(filter_parameters);
            app = app.compute_active_conductances();
            app = app.compute_active_currents();
            app = app.compute_leakage_currents();
            app = app.compute_membrane_currents();
            app = app.zero_phase_filter_Im(filter_parameters);
            app = app.compute_passive_conductances_nopage();
            app = app.compute_passive_currents();
            app = app.compute_stats();
       end

       function app = dynamics_plots(app)
            tStart = tic;
            app(1).fig = figure('Name', strcat(app(1).filename, ' Dynamics plots'));
            tiledlayout(4, length(app));
            ax = cell(4, length(app));
            for m = 1: 1: 4
                for k = 1: 1: length(app)
                    ax{m, k} = nexttile;
                    if m == 1
                        plot(app(k).Vm, app(k).Im);
                        xlabel('Vm (V)');
                        if k == 1
                            ylabel('Im (A)');
                        end
                        title(strcat(app(k).condition'));
                    elseif m == 2
                        plot(app(k).Vm, app(k).Iactive);
                        xlabel('Vm (V)');
                        if k == 1
                            ylabel('Iactive (A)');
                        end
                    elseif m == 3
                        plot(app(k).Vm, app(k).Ileak);
                        xlabel('Vm (V)');
                        if k == 1
                            ylabel('Ileak (A)');
                        end
                    elseif m == 4
                        plot(app(k).Im, app(k).Iactive);
                        xlabel('Im (A)');
                        if k == 1
                            ylabel('Iactive (A)');
                        end
                    end
                end
            end
            fprintf('[%d secs] Dynamics Plotting data\n', toc(tStart));
       end
   end
end