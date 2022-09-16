function [meta_stats, stats] = analyze(filename, paradigms, response_times)
    %% Data reading parameters.
    rec = WholeCellRecording2(filename, paradigms, response_times);
    %% Filter parameters.
    filter_parameters.CutOffFrequency = 150; % Vm cuttoff. The lower this value, the smoother the traces get, but stay above 80 (Hz) for now.
    filter_parameters.CutOffFrequency2 = 20; % Im cuttoff.
    filter_parameters.FilterOrder = 100;
    filter_parameters.PassbandRipple = 0.01;
    filter_parameters.StopbandAttenuation = 80;
    %% Analysis.
    rec = rec.call(filter_parameters);
    %% Plotting.
    rec = rec.plots();
    %% Computing stats.
    stats = rec.get_stats();
    %% Computing meta-stats.
    meta_stats = rec.compute_meta_stats();
end