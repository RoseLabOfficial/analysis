clc; clear;
filename = "./example_data.xlsx";
paradigms = ["1pulse", "80pps-2pulses", "80pps-4pulses", "80pps-6pulses", "80pps-8pulses"];
response_durations = [0.6, 0.6, 0.6, 0.6, 0.6];

membrane_potential_filters.firdifferentiator.Fp = 80;
membrane_potential_filters.firdifferentiator.Fst = 100;
membrane_potential_filters.firdifferentiator.Ap = 0.01;
membrane_potential_filters.firdifferentiator.Ast = 65;

membrane_potential_filters.firlowpass.Fp = 400;
membrane_potential_filters.firlowpass.Fst = 500;
membrane_potential_filters.firlowpass.Ap = 0.01;
membrane_potential_filters.firlowpass.Ast = 65;

ds = WholeCellRecordingV2(filename, paradigms, response_durations, membrane_potential_filters);
ds = ds.call();
ds = ds.plot_estimations();