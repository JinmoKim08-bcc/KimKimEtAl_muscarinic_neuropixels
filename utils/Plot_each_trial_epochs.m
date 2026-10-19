function Plot_each_trial_epochs(trials, verbose)

if verbose == 1 % do plot
    num_trials = size(trials, 2);
    lickportsIN = ones(size(trials, 2), 1);
    lickON = ones(size(trials, 2), 1);

    for trial_idx = 1:size(trials, 2)
        start_time = trials(trial_idx).time(1);
        portsIN_idx = find(strcmp(trials(trial_idx).logs(:, 2), 'lickportsIN') == 1);
        lickportsIN(trial_idx, 1) = str2double(trials(trial_idx).logs(portsIN_idx, 1)) - start_time;
        lickON(trial_idx, 1) = str2double(trials(trial_idx).logs(portsIN_idx+1, 1)) - start_time;
    end

    % figure 
    figure;
    hold on;

    % raster plot 
    for trial = 1:num_trials
        % lickportsIN 
        x = lickportsIN(trial);
        y = trial;

        plot([x, x], [y-0.4, y+0.4], 'b', 'LineWidth', 1);

        x2 = lickON(trial);
        plot([x2, x2], [y-0.4, y+0.4], 'r', 'LineWidth', 1);
    end

    xlabel('Time (s)');
    ylabel('Trial');
    ylim([0, num_trials+1]);
    xlim([min([lickportsIN; lickON]) - 0.1, max([lickportsIN; lickON]) + 0.1]);
    xline(0.05, 'g:', 'LineWidth', 1.5)
    xline(2.05, 'g:', 'LineWidth', 1.5)

    legend({'lickportsIN','lickOn'}, 'Location', 'best');

    title('Raster plot of lickportsIN and lickOn events');
    set(gca, 'FontSize', 15)
    hold off;
end


end

