function plot_reverse_estimation(ds, exe, inh)
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
                        plot(exe(i, j).time, exe(i, j).Vm);
                        ax{j, k}.Title.String = exe(i, j).paradigm;
                        if j == 1
                            ax{j, k}.YLabel.String = 'Vm (V)';
                        end
                    case 2
                        plot(exe(i, j).time, exe(i, j).ge, 'r', exe(i, j).time, exe(i, j).gi, 'b');
                        if j == 1
                            ax{j, k}.YLabel.String = 'G (S)';
                        end
    %                     ax{j, k}.XLabel.String = 'time (sec)'; 
                    case 3
                        plot(inh(i, j).time, inh(i, j).Vm);
                        ax{j, k}.Title.String = inh(i, j).paradigm;
                        if j == 1
                            ax{j, k}.YLabel.String = 'Vm (V)';
                        end
                    case 4
                        plot(inh(i, j).time, inh(i, j).ge, 'r', inh(i, j).time, inh(i, j).gi, 'b');
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
end