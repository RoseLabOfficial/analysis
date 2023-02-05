%% Data reading parameters.
filename = "example_data.xlsx"; % The path+name of the file goes here.
paradigms = ["10pps", "20pps", "40pps", "60pps", "80pps"]; % Different rates in the file go here; these should be the same as the worksheets in the file.
active_duration = [0.6, 0.4, 0.35, 0.3, 0.3]; % How long does the response to simulus last at each pulse rate? (secs).
%% Read
n = NeuronObj(filename, paradigms, active_duration);