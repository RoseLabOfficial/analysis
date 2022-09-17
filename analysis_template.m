%% Data reading parameters.
filename = "example_data.xlsx"; % The path+name of the file goes here.
paradigms = ["10pps", "20pps", "40pps", "60pps", "80pps"]; % Different rates in the file go here; these should be the same as the worksheets in the file.
active_duration = [0.6, 0.4, 0.35, 0.3, 0.3]; % How long does the response to simulus last at each pulse rate? (secs).
%% Run analysis
[ds, stats, meta_stats] = analyze(filename, paradigms, active_duration, [], 0); %[] indicates default filter settings. 1 = overwrite stats or 0 = create new stats sheet.
% ds: data structure, WholeCellRecording object.
% stats: stats variable with stats for each paradigm.
% meta_stats: stats variable with meta stats over all paradigms.