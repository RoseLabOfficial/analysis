%% Import data
filename = "./example_data.xlsx";
paradigms = ["1pulse", "80pps-2pulses", "80pps-4pulses", "80pps-6pulses", "80pps-8pulses"];
response_durations = [0.6, 0.6, 0.6, 0.6, 0.6];
%% Filters
filters = {};
filters.lowpass = {};
filters.lowpass.Fp = 100;
filters.lowpass.Fst = 140;
filters.lowpass.Ap = 0.01;
filters.lowpass.Ast = 60;
filters.notch = {};
filters.notch.Fst = 60;
filters.notch.BW = 6;
%% Whole-cell recordings data structure
ds = WholeCellRecordingV2(filename, paradigms, response_durations, filters);
%% Analyze
ds = ds.call();
%% Compute Stats
ds = ds.compute_means(1);
meta_stats = ds.compute_meta_stats();
disp(meta_stats);
%% Build Plots
ds = ds.build_plots();
% ds = ds.plot_estimations();
ds = ds.plot_meta_stats(meta_stats);