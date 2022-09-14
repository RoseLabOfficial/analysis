filename = "example_data.xlsx"; % The path+name of the file goes here.
rates = ["10pps", "20pps", "40pps"]; % Different rates in the file go here; these should be the same as the worksheets in the file.
active_duration = [0.6, 0.4, 0.4]; % How long does the response to simulus last at each pulse rate? (secs).
filter_parameters.CutOffFrequency = 150; % The lower this value, the smoother the traces get, but stay above 80 (Hz) for now.
filter_parameters.FilterOrder = 100;
filter_parameters.PassbandRipple = 0.01;
filter_parameters.StopbandAttenuation = 80;
data = WholeCellRecording(filename, rates, active_duration);
data = data.zero_phase_filter_Vm(filter_parameters);
data = data.compute_active_conductances();
data = data.compute_active_currents();
data = data.compute_leakage_currents();
data = data.compute_membrane_currents();
filter_parameters.CutOffFrequency = 20;
data = data.zero_phase_filter_Im(filter_parameters);
data = data.compute_passive_conductances();
% rec = rec.plots();
[data, stats] = data.generate_stats();
fprintf("Stats for: %s \n", filename); 
disp(stats);
% rec = rec.dynamics_plots();
figure('Name', 'WholeCellRecording');
plot(data(1).time, data(1).ge, 'r', data(1).time, data(1).gi, 'b');