classdef WholeCellRecording
   properties
       %% Sampling/Acquistion properties.
       Fs % Sampling Frequency of acquisition (Hz).
       response_samples % How long the response lasts in a sampling domain (samples).
       
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
       time % Pseudo time array created as per Fs, to analyze the temporal dynamics of excitation and inhibition.
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
       condition % condition (string) is different types of stimuli, for pulse rate 5pps, 10pps, 60pps, etc., or for duration 20ms, 40ms, 160ms, etc.
       filename % filename (string) is the name of the file containing the current clamp data along with various parameters for cell and membrane potential constants.
       
   end
   properties(Access=private)
       %% fig properties
       fig
   end
   %% Methods
   methods
       %% Constructor
       function DATA = WholeCellRecording(filename, conditions, response_durations)
            if nargin > 0
                tStart = tic;
                %% Input format check
                if ~isa(filename, 'string')
                    error('filename must be a string'); 
                end
                if ~isa(conditions, 'string')
                    error('conditions (rates/durations) must be strings'); 
                end
                if ~isa(response_durations, 'double')
                    error('responseDurations must be doubles');
                end
                if ~isequal(size(conditions), size(response_durations))
                    error('Dimensions of conditions and responseDurations must match.');
                end
                DATA(1).filename = filename;
                %% Reading worksheets from excel files.
                for k = numel(conditions):-1:1
                    data = xlsread(filename, conditions(k));
                    [samples, ~] = size(data);
                    DATA(k).condition = conditions(k);
                    DATA(k).time = data(2:samples, 1);
                    DATA(k).Fs = 1/(DATA(k).time(2, 1) - DATA(k).time(1, 1));
                    DATA(k).response_samples = floor(DATA(k).Fs*response_durations(k));
                    if DATA(k).response_samples > samples-1
                       DATA(k).response_samples = samples-1;
                    end
                    DATA(k).Vm = data(2:samples, 2:end)*1e-3;
                    DATA(k).Iactive = zeros(size(DATA(k).Vm));
                    DATA(k).ge = zeros(size(DATA(k).Vm, 1), 1);
                    DATA(k).gi = zeros(size(DATA(k).Vm, 1), 1);
                    DATA(k).Ie = zeros(size(DATA(k).Vm));
                    DATA(k).Ii = zeros(size(DATA(k).Vm));
                    parameters = xlsread(filename, "parameters_"+conditions(k));
                    DATA(k).Iinj = parameters(1:end, 1);
                    DATA(k).Cm = parameters(1:end, 2);
                    DATA(k).Rin = parameters(1:end, 3);
                    DATA(k).Er = parameters(1:end, 4);
                    DATA(k).Ee = parameters(1:end, 5);
                    DATA(k).Ei = parameters(1:end, 6);
                    DATA(k).Et = parameters(1:end, 7);
                    DATA(k).Eact = parameters(1:end, 8);
                    DATA(k).Ess = parameters(1:end, 9);
                    DATA(k).xalpha = parameters(1:end, 10);
                    DATA(k).xbeta = parameters(1:end, 11);
                    try
                        DATA(k).sps = parameters(1:end, 12);
                    catch
                        DATA(k).sps = zeros(size(parameters, 1), 1);
                    end
                    try
                        DATA(k).Eref = parameters(1:end, 13);
                    catch
                        DATA(k).Eref = DATA(k).Ess;
                    end
                    for m = 1: 1: size(DATA(k).Iinj, 1)
                        DATA(k).Iinj(m) = (1/DATA(k).Rin(m)).*(DATA(k).Ess(m) - DATA(k).Er(m));
                        DATA(k).Vm(:, m) = DATA(k).Vm(:, m) - DATA(k).Vm(1, m) + DATA(k).Ess(m);
                    end
                end
                fprintf('[%d secs] Read %s\n', toc(tStart), filename);
            end
       end
       function DATA = zero_phase_filter_Vm(DATA, filter_parameters)
           tStart = tic;
           Fnorm = filter_parameters.CutOffFrequency/(DATA(1).Fs/2);
           lpFilt = designfilt('lowpassiir', ...
                                'PassbandFrequency', Fnorm, ...
                                'FilterOrder', filter_parameters.FilterOrder, ...
                                'PassbandRipple', filter_parameters.PassbandRipple, ...
                                'StopbandAttenuation', filter_parameters.StopbandAttenuation);
           for k = 1: 1: length(DATA)
               [~, clamps] = size(DATA(k).Vm);
               v = cat(1, ones(size(DATA(k).Vm)).*(DATA(k).Vm(1, :)), DATA(k).Vm);
               v = cat(1, v, ones(size(DATA(k).Vm)).*(DATA(k).Vm(end, :)));
               vn = v - v(1, :);
               Vmf = zeros(size(vn));
               for clamp = 1: 1: clamps
                   Vmf(:, clamp) = filtfilt(lpFilt, vn(:, clamp));
               end
               DATA(k).Vm = Vmf(size(DATA(k).Vm, 1)+1:end-size(DATA(k).Vm, 1), :) + DATA(k).Vm(1, :);
           end
           fprintf('[%d secs] Zero phase filtering Vm \n', toc(tStart));
       end
       function DATA = compute_active_conductances(DATA)
           tStart = tic;
           for k = 1: 1: length(DATA)
               DATA(k).alpha = ((1./DATA(k).Rin)./(2.*(DATA(k).Eact-DATA(k).Ess)));
               DATA(k).beta = DATA(k).alpha.*(DATA(k).Et - DATA(k).Ess);
               DATA(k).alpha = DATA(k).alpha.*DATA(k).xalpha;
               DATA(k).beta = DATA(k).beta.*DATA(k).xbeta;
           end
           fprintf('[%d secs] Computed active conductances (constants)\n', toc(tStart));
       end
       function DATA = compute_active_currents(DATA)
           tStart = tic;
           for k = 1: 1: length(DATA)
               DATA(k).Iactive = zeros(size(DATA(k).Vm));
              for m = 1: 1: size(DATA(k).Iactive, 2)
                  DATA(k).Iactive(:, m) = (DATA(k).alpha(m).*(DATA(k).Vm(:, m)-DATA(k).Et(m)).*(DATA(k).Vm(:, m)-DATA(k).Er(m))) ...
                      + (DATA(k).beta(m).*(DATA(k).Vm(:, m)-DATA(k).Er(m)));
              end
           end
           fprintf('[%d secs] Computed active currents\n', toc(tStart));
       end
       function DATA = compute_leakage_currents(DATA)
           tStart = tic;
           for k = 1: 1: length(DATA)
               DATA(k).Ileak = zeros(size(DATA(k).Vm));
              for m = 1: 1: size(DATA(k).Iactive, 2)
                  DATA(k).Ileak(:, m) = (1./DATA(k).Rin(m)).*(DATA(k).Vm(:, m)-DATA(k).Er(m));
              end
           end
           fprintf('[%d secs] Computed leakage currents\n', toc(tStart));
       end
       function DATA = compute_membrane_currents(DATA)
           tStart=tic;
           for k = 1:1:length(DATA)
               DATA(k).Im = zeros(size(DATA(k).Vm));
               for m = 1: 1: size(DATA(k).Vm, 2)
                   Vm_appended = cat(1, DATA(k).Vm(1, m), DATA(k).Vm(:, m), DATA(k).Vm(end, m));
                   dVmdt = diff(Vm_appended);
                   DATA(k).Im(:, m) = DATA(k).Cm(m).*dVmdt(1:end-1).*(DATA(k).Fs);
               end
           end
           fprintf('[%d secs] Computed membrane currents\n', toc(tStart));
       end
       function DATA = zero_phase_filter_Im(DATA, filter_parameters)
            tStart = tic;
            Fnorm = filter_parameters.CutOffFrequency/(DATA(1).Fs/2);
            lpFilt = designfilt('lowpassiir', ...
                                'PassbandFrequency', Fnorm, ...
                                'FilterOrder', filter_parameters.FilterOrder, ...
                                'PassbandRipple', filter_parameters.PassbandRipple, ...
                                'StopbandAttenuation', filter_parameters.StopbandAttenuation);
            for k = 1: 1: length(DATA)
                [~, clamps] = size(DATA(k).Im);
                i = cat(1, ones(size(DATA(k).Im)).*(DATA(k).Im(1, :)), DATA(k).Im);
                i = cat(1, i, ones(size(DATA(k).Im)).*(DATA(k).Im(end, :)));
                in = i - i(1, :);
                Imf = zeros(size(in));
                for clamp = 1: 1: clamps
                    Imf(:, clamp) = filtfilt(lpFilt, in(:, clamp));
                end
                DATA(k).Im = Imf(size(DATA(k).Im, 1)+1:end-size(DATA(k).Im, 1), :) + DATA(k).Im(1, :);
            end
            fprintf('[%d secs] Zero phase filtering Im \n', toc(tStart));
       end
       function DATA = compute_passive_conductances(DATA)
           tStart = tic;
           for k = 1: 1: length(DATA)
               for n = 1: 1: size(DATA(k).time, 1)
                   A = zeros(2, 2);
                   B = zeros(2, 1);
                   A(1, 1) = sum((DATA(k).Vm(n, :) - DATA(k).Ee').^2);
                   A(1, 2) = sum((DATA(k).Vm(n, :) - DATA(k).Ee').*(DATA(k).Vm(n, :) - DATA(k).Ei'));
                   A(2, 1) = A(1, 2);
                   A(2, 2) = sum((DATA(k).Vm(n, :) - DATA(k).Ei').^2);
                   C = DATA(k).Im(n, :) - DATA(k).Iinj' - DATA(k).Iactive(n, :) + DATA(k).Ileak(n, :);
                   B(1, 1) = -sum(C.*(DATA(k).Vm(n, :) - DATA(k).Ee'));
                   B(2, 1) = -sum(C.*(DATA(k).Vm(n, :) - DATA(k).Ei'));
                   G = pinv(A)*B;
                   DATA(k).ge(n, 1) = G(1, 1);
                   DATA(k).gi(n, 1) = G(2, 1);
               end
               DATA(k).Ie = DATA(k).ge.*(DATA(k).Vm - DATA(k).Ee');
               DATA(k).Ii = DATA(k).gi.*(DATA(k).Vm - DATA(k).Ei');
           end
           fprintf('[%d secs] Computed passive conductances\n', toc(tStart));
       end
       function DATA = plots(DATA)
            tStart = tic;
            DATA(1).fig = figure('Name', strcat(DATA(1).filename, ' Reconstructions'));
            tiledlayout(6, length(DATA));
            ax = cell(6, length(DATA));
            for m = 1: 1: 6
               for k = 1: 1: length(DATA)
                   ax{m, k} = nexttile;
                   if m == 1
                      plot(DATA(k).time, DATA(k).Vm);
                      if k == 1
                        ylabel('Vm (V)');
                      end
                      title(strcat(DATA(k).condition));
                   elseif m == 2
                      plot(DATA(k).time, DATA(k).Iactive);
                      if k == 1
                        ylabel('Iactive (A)');
                      end
                   elseif m == 3
                      plot(DATA(k).time, DATA(k).Ileak);
                      if k == 1
                        ylabel('Ileak (A)');
                      end
                   elseif m == 4
                      plot(DATA(k).time, DATA(k).Im);
                      if k == 1
                        ylabel('Im (A)');
                      end
                   elseif m == 5
                       plot(DATA(k).time, DATA(k).ge, 'r', DATA(k).time, DATA(k).gi, 'b');
                       if k == 1
                         ylabel('G (S)');
                       end
                       xlabel('Time (sec)');
                   elseif m == 6
                       plot(DATA(k).time, -1.*DATA(k).Ie(:, 1), 'r', DATA(k).time, DATA(k).Ii(:, 1), 'b');
                       if k == 1
                           ylabel('I (A)');
                       end
                       xlabel('Time (sec)');
                   end
               end
               linkaxes([ax{m, :}], 'xy');
            end
            fprintf('[%d secs] Plotting data\n', toc(tStart));
       end
       function DATA = dynamics_plots(DATA)
            tStart = tic;
            DATA(1).fig = figure('Name', strcat(DATA(1).filename, ' Dynamics plots'));
            tiledlayout(4, length(DATA));
            ax = cell(4, length(DATA));
            for m = 1: 1: 4
                for k = 1: 1: length(DATA)
                    ax{m, k} = nexttile;
                    if m == 1
                        plot(DATA(k).Vm, DATA(k).Im);
                        xlabel('Vm (V)');
                        if k == 1
                            ylabel('Im (A)');
                        end
                        title(strcat(DATA(k).condition'));
                    elseif m == 2
                        plot(DATA(k).Vm, DATA(k).Iactive);
                        xlabel('Vm (V)');
                        if k == 1
                            ylabel('Iactive (A)');
                        end
                    elseif m == 3
                        plot(DATA(k).Vm, DATA(k).Ileak);
                        xlabel('Vm (V)');
                        if k == 1
                            ylabel('Ileak (A)');
                        end
                    elseif m == 4
                        plot(DATA(k).Im, DATA(k).Iactive);
                        xlabel('Im (A)');
                        if k == 1
                            ylabel('Iactive (A)');
                        end
                    end
                end
            end
            fprintf('[%d secs] Dynamics Plotting data\n', toc(tStart));
       end
       function [DATA, stats] = generate_stats(DATA)
           tStart = tic;
           net_conductances = zeros(length(DATA), 2);
           mean_conductances = zeros(length(DATA), 2);
           mean_polarizations = zeros(length(DATA), 2);
           mean_active_current = zeros(length(DATA), 1);
           spikes_per_rep = zeros(length(DATA), 1);
           pulse_rates = strings([length(DATA), 1]);
           data_count = zeros(length(DATA), 1);
           DATA(1).fig = figure('Name', strcat(DATA(1).filename, ' Stats'));
           set(DATA(1).fig,'defaultAxesColorOrder',[[0.5 0.5 0]; [0 0 0]]);
           tiledlayout(2, 1);
           ax = cell(2, 1);
           for k = 1: 1: length(DATA)
               del_Vm = (DATA(k).Vm(:, 1) - DATA(k).Eref(1));
               DATA(k).depolarizations = del_Vm(:, 1).*(del_Vm(:, 1)>0);
               DATA(k).hyperpolarizations = del_Vm(:, 1).*(del_Vm(:, 1)<0);
               DATA(k).excitation = DATA(k).ge.*(DATA(k).ge>0);
               DATA(k).inhibition = DATA(k).gi.*(DATA(k).gi>0);
               resultant_conductance = DATA(k).excitation(1:DATA(k).response_samples) - DATA(k).inhibition(1:DATA(k).response_samples);
               net_conductances(k, 1) = mean(resultant_conductance.*(resultant_conductance>0));
               net_conductances(k, 2) = mean(resultant_conductance.*(resultant_conductance<0));
               mean_conductances(k, 1) = mean(DATA(k).excitation(1:DATA(k).response_samples));
               mean_conductances(k, 2) = -1.*mean(DATA(k).inhibition(1:DATA(k).response_samples));
               mean_polarizations(k, 1) = mean(DATA(k).depolarizations(1:DATA(k).response_samples));
               mean_polarizations(k, 2) = mean(DATA(k).hyperpolarizations(1:DATA(k).response_samples));
               mean_active_current(k, 1) = mean(DATA(k).Iactive(1:DATA(k).response_samples, 1));
               spikes_per_rep(k, 1) = DATA(k).sps(1);
               pulse_rates(k, 1) = DATA(k).condition;
               data_count(k, 1) = k;
%                fprintf('mean depolarization %d, mean hyperpolarization %d, mean excitation %d, mean inhibition %d\n', mean(DATA(k).depolarizations), mean(DATA(k).hyperpolarizations), mean(DATA(k).excitation), mean(DATA(k).inhibition));
           end
           ax{1, 1} = nexttile;
           yyaxis left;
           mean_plt = bar(data_count, mean_conductances, 0.3, 'stacked', 'FaceColor', 'flat');
           hold on;
           net_plt = bar(data_count, net_conductances, 0.5, 'stacked', 'FaceColor', 'flat');
           hold off;
           ylabel('Mean Conductance (S)');
           xlabel('Pulse Rate');
           for k = 1: 1: length(DATA)
               mean_plt(1).CData(k, :) = [1, 0, 0];
               mean_plt(2).CData(k, :) = [0, 0, 1];
               net_plt(1).CData(k, :) = [1, 0, 0];
               net_plt(2).CData(k, :) = [0, 0, 1];
           end
           y1max = max(mean_conductances, [], 'all');
           y1max = y1max+0.1*y1max;
           if y1max == 0
               ylim([0, 1]);
           else
               ylim([-y1max, y1max]);
           end
           y2max = max(spikes_per_rep, [], 'all');
           y2max = y2max+0.1*y2max;
           yyaxis right;
           plot(data_count, spikes_per_rep, '-ok');
           ylabel('Spikes per stim. rep. (SPS)');
           if y2max == 0
               ylim([0, 1]);
           else
               ylim([-y2max, y2max]);
           end
           set(ax{1, 1},'xticklabel',pulse_rates);
           ax{2, 1} = nexttile;
           yyaxis left;
           mean_plt = bar(data_count, mean_polarizations, 'stacked', 'FaceColor', 'flat');
           ylabel('Mean Polarizaions (V)');
           xlabel('Pulse Rate');
           for k = 1: 1: length(DATA)
               mean_plt(1).CData(k, :) = [1, 0, 0];
               mean_plt(2).CData(k, :) = [0, 0, 1];
           end
           y1max = max(mean_polarizations, [], 'all');
           y1max = y1max+0.1*y1max;
           if y1max == 0
               ylim([0, 1]);
           else
               ylim([-y1max, y1max]);
           end
           y2max = max(spikes_per_rep, [], 'all');
           y2max = y2max+0.1*y2max;
           yyaxis right;
           plot(data_count, spikes_per_rep, '-ok');
           ylabel('Spikes per stim. rep. (SPS)');
           if y2max == 0
               ylim([0, 1]);
           else
               ylim([-y2max, y2max]);
           end
           set(ax{2, 1},'xticklabel',pulse_rates);
           fprintf('[%d secs] Plotting stats\n', toc(tStart));
           stats = table(pulse_rates, spikes_per_rep, mean_polarizations, mean_conductances, net_conductances, mean_active_current);
       end
   end
end