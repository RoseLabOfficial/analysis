function [ds, stats, meta_stats] = analyze(filename, paradigms, response_times, filter_parameters, overwrite_stats)
    %% Data reading parameters.
    ds = WholeCellRecording(filename, paradigms, response_times);
    %% Filter parameters.
    if nargin < 4 || isempty(filter_parameters)
        filter_parameters.CutOffFrequency = 150; % Vm cuttoff. The lower this value, the smoother the traces get, but stay above 80 (Hz) for now.
        filter_parameters.CutOffFrequency2 = 20; % Im cuttoff.
        filter_parameters.FilterOrder = 100;
        filter_parameters.PassbandRipple = 0.01;
        filter_parameters.StopbandAttenuation = 80;
    end
    if nargin < 5 || isempty(overwrite_stats)
        overwrite_stats = 0;
    end
    %% Analysis.
    ds = ds.call(filter_parameters);
    %% Computing stats.
    stats = ds.get_stats();
    disp(stats);
    %% Computing meta-stats.
    meta_stats = ds.get_meta_stats();
    disp(meta_stats);
    %% Writing stats to file.
    if overwrite_stats == 0
        ds.write_stats_to_file(stats, filename, strcat('Stats', strrep(datestr(now), ':', '-')));
    else
        ds.write_stats_to_file(stats, filename, 'Stats');
    end
    %% Writing meta-stats to file.
    if overwrite_stats == 0
        ds.write_meta_stats_to_file(meta_stats, filename, strcat('MetaStats', strrep(datestr(now), ':', '-')));
    else
        ds.write_meta_stats_to_file(meta_stats, filename, 'MetaStats');
    end
    %% Plotting.
    ds = ds.plots();
end