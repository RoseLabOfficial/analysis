%% Data reading parameters.
filename = "example_data.xlsx"; % The path+name of the file goes here.
paradigms = ["10pps", "20pps", "40pps", "60pps", "80pps"]; % Different rates in the file go here; these should be the same as the worksheets in the file.
active_duration = [0.6, 0.4, 0.35, 0.3, 0.3]; % How long does the response to simulus last at each pulse rate? (secs).
%% Run analysis
[ds, stats, meta_stats] = analyze(filename, paradigms, active_duration, [], 0); %[] indicates default filter settings. 1 = overwrite stats or 0 = create new stats sheet.
% ds: data structure, WholeCellRecording object.
% stats: stats variable with stats for each paradigm.
% meta_stats: stats variable with meta stats over all paradigms.
return
%% Reverse Estimate Vm
ds_exe = ds;
ds_exe = ds_exe.reverse_estimate_membrane_potential(1, 0);
ds_inh = ds;
ds_inh = ds_inh.reverse_estimate_membrane_potential(0, 1);
%% Plot reverse estimations
[m, n] = size(ds);
nplots = 4;
figure('Name', strcat(ds(1, 1).filename, 'Reverse Estimations'));
for i = 1: 1: m
    tiledlayout(nplots, n);
    ax = cell(nplots, n);
    for k = 1: 1: nplots
        for j = 1: 1: n
            ax{j, k} = nexttile;
            switch k
                case 1
                    plot(ds_exe(i, j).time, ds_exe(i, j).Vm);
                    ax{j, k}.Title.String = ds_exe(i, j).paradigm;
                    if j == 1
                        ax{j, k}.YLabel.String = 'Vm (V)';
                    end
                case 2
                    plot(ds_exe(i, j).time, ds_exe(i, j).ge, 'r', ds_exe(i, j).time, ds_exe(i, j).gi, 'b');
                    if j == 1
                        ax{j, k}.YLabel.String = 'G (S)';
                    end
%                     ax{j, k}.XLabel.String = 'time (sec)'; 
                case 3
                    plot(ds_inh(i, j).time, ds_inh(i, j).Vm);
                    ax{j, k}.Title.String = ds_inh(i, j).paradigm;
                    if j == 1
                        ax{j, k}.YLabel.String = 'Vm (V)';
                    end
                case 4
                    plot(ds_inh(i, j).time, ds_inh(i, j).ge, 'r', ds_inh(i, j).time, ds_inh(i, j).gi, 'b');
                    if j == 1
                        ax{j, k}.YLabel.String = 'G (S)';
                    end
                    ax{j, k}.XLabel.String = 'time (sec)'; 
            end
            linkaxes([ax{j, :}], 'x');
        end
        linkaxes([ax{:, k}], 'y');
    end
end
