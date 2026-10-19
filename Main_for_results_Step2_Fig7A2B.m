clc
clear
close all

%% Hyperparameters to adjust
npy_path = "utils\npy-matlab\npy-matlab\";
addpath(npy_path)

param = struct();

save_path = 'C:\Users\sirom\Documents\KimKimEtAl_muscarnic_neuropixels\';
parts_path = '\fig6\';

if ~isfolder(strcat(save_path, parts_path))
    mkdir(strcat(save_path, parts_path));
end

%% Hyperparameters that most do not need to be adjusted
param.bin_time = 0.001; % [sec]; temporal binning size to calculate firing rate
param.ITI_time = 2; % [sec]; Inter-trial interval duration
param.sample_time = 0.1; % [sec]; timing given a sample
param.delay_time = 2; % [sec]
param.response_time = 0.5; % [sec]
param.t_max = param.sample_time + param.delay_time + param.response_time; % [sec]; maximal end time to analyze

% parameters for direction selectivity (RR vs LL)
param.epoch_names = ["Stimulus"; "Delay Early"; "Delay Late"; "Response"];
param.epochs = [0, 0.1; % Stimulus [sec]
    0.1, 1.1; % Delay Early
    1.1, 2.1; % Delay Late
    2.1, 2.6]; % Response
param.delay_epoch_idx = [2, 3]; % Delay Early, Delay Late
param.smoothing_pts = 100; % [idx]; trailing moving average
param.sel_alpha = 0.05; % raw p-value of ranksum (no correction)
param.num_sampled_units = 100; % units sampled with replacement per repeat

% parameters for LDA with temporal resolution
param.win_time = 0.2; % [sec]; causal window
param.step_time = 0.05; % [sec]
param.start_time = 1; % [sec]; first window starts 1 sec after the ITI onset
param.num_repeat = 100;
param.num_chosen_units = 100;
param.num_replaced_list = [0, 10, 20, 10, 20]; % Atropine units replaced by Before units
param.pool_list = [0, 1, 1, 2, 2]; % 1: all Before units; 2: Before delay late selective units
param.num_trials_per_class = 50; % correct trials per class
param.num_folds = 5;
param.rng_repeat = 10000; % offset of rng seeds
param.t_mark = [0.8, 1.5]; % [sec]; ROC curves are kept at these times
param.bar_epoch = [1.1, 2.1]; % [sec]; delay late window for the auROC bar plot
param.mix_pairs = [1, 2; 1, 3; 1, 4; 1, 5; 2, 4; 3, 5]; % compared populations in the bar plot

param.save_path = save_path;
param.version = "v16fd";

cond_name_list = ["Before"; "Atropine"];
num_sess_list = [5, 4];
trial_types = ["RR"; "LL"; "LR"; "RL"];

mix_labels = ["A100";
    "A90 + Before random 10 units";
    "A80 + Before random 20 units";
    "A90 + Before delay late direction selective random 10 units";
    "A80 + Before delay late direction selective random 20 units"];

time_bin = -param.ITI_time:param.bin_time:param.t_max-param.bin_time;
num_pts = numel(time_bin);
ITI_dur_bin = round(param.ITI_time/param.bin_time);

% causal windows for LDA
win_pts = round(param.win_time/param.bin_time);
step_pts = round(param.step_time/param.bin_time);
ub_list = (round(param.start_time/param.bin_time)+win_pts):step_pts:num_pts;
lb_list = ub_list-win_pts+1;
num_win = numel(ub_list);
time_axis = time_bin(ub_list)*1000; % [ms]; end of each causal window

t_mark_idx = zeros(size(param.t_mark));
for mark_idx = 1:numel(param.t_mark)
    [~, t_mark_idx(mark_idx)] = min(abs(time_axis-param.t_mark(mark_idx)*1000));
end

%% Read calculated firing vectors
load(strcat(param.save_path, 'results_', param.version, '.mat'), 'results');

%% Direction selectivity of each unit
% ITI subtraction per trial -> trailing 100-ms mean -> epoch mean -> RR vs LL ranksum
selectivity = struct();

for cond_idx = 1:numel(cond_name_list)
    p_raw = [];
    auroc = [];

    for sess_idx = 1:num_sess_list(cond_idx)
        fprintf('Cond: [%d/%d], Sess: [%d/%d]\n', cond_idx, numel(cond_name_list), sess_idx, num_sess_list(cond_idx));

        vals_R = results(cond_idx).sess(sess_idx).RR.allTrials.vals; % units x time x trials
        vals_R = vals_R(:, ITI_dur_bin+1:num_pts, :) - mean(vals_R(:, 1:ITI_dur_bin, :), 2);
        vals_R = movmean(vals_R, [param.smoothing_pts-1, 0], 2, 'Endpoints', 'shrink');

        vals_L = results(cond_idx).sess(sess_idx).LL.allTrials.vals;
        vals_L = vals_L(:, ITI_dur_bin+1:num_pts, :) - mean(vals_L(:, 1:ITI_dur_bin, :), 2);
        vals_L = movmean(vals_L, [param.smoothing_pts-1, 0], 2, 'Endpoints', 'shrink');

        num_units = size(vals_R, 1);
        temp_p = ones(num_units, size(param.epochs, 1));
        temp_auroc = zeros(num_units, size(param.epochs, 1));

        for epoch_idx = 1:size(param.epochs, 1)
            epoch_range = round(param.epochs(epoch_idx, 1)/param.bin_time)+1:round(param.epochs(epoch_idx, 2)/param.bin_time);
            x_R = reshape(mean(vals_R(:, epoch_range, :), 2), num_units, []); % units x trials
            x_L = reshape(mean(vals_L(:, epoch_range, :), 2), num_units, []);

            for unit_idx = 1:num_units
                x = x_R(unit_idx, :);
                y = x_L(unit_idx, :);

                if all([x, y] == x(1))
                    temp_auroc(unit_idx, epoch_idx) = 0.5; % identical values: p = 1
                else
                    temp_p(unit_idx, epoch_idx) = ranksum(x, y);
                    ranks = tiedrank(double([x, y]'));
                    temp_auroc(unit_idx, epoch_idx) = (sum(ranks(1:numel(x))) - numel(x)*(numel(x)+1)/2)/(numel(x)*numel(y));
                end
            end
        end

        p_raw = [p_raw; temp_p];
        auroc = [auroc; temp_auroc];
    end

    selectivity(cond_idx).p_raw = p_raw;
    selectivity(cond_idx).auroc = auroc;
    selectivity(cond_idx).is_selective = p_raw < param.sel_alpha;
end

for cond_idx = 1:numel(cond_name_list)
    fprintf('%s: %d units, selective units (Stimulus, Delay Early, Delay Late, Response) = %d, %d, %d, %d\n', ...
        cond_name_list(cond_idx), size(selectivity(cond_idx).is_selective, 1), sum(selectivity(cond_idx).is_selective, 1));
end

%% Causal window features for LDA
% ITI-subtracted mean firing rate in each 200-ms causal window (units x windows x trials)
% Atropine sessions first, then Before sessions
features = struct();
num_units_list = [];
num_incor_trials = inf;

temp_sess_idx = 0;
for cond_idx = [2, 1]
    for sess_idx = 1:num_sess_list(cond_idx)
        temp_sess_idx = temp_sess_idx + 1;

        for trial_idx = 1:numel(trial_types)
            vals = results(cond_idx).sess(sess_idx).(trial_types(trial_idx)).allTrials.vals;
            vals_baseline = mean(vals(:, 1:ITI_dur_bin, :), 2);
            vals_cumsum = cumsum(vals, 2);

            temp_features = zeros(size(vals, 1), num_win, size(vals, 3), 'single');
            for win_idx = 1:num_win
                temp_features(:, win_idx, :) = (vals_cumsum(:, ub_list(win_idx), :) - vals_cumsum(:, lb_list(win_idx)-1, :))/win_pts - vals_baseline;
            end
            features(temp_sess_idx).(trial_types(trial_idx)) = temp_features;
        end

        num_units_list = [num_units_list; size(vals, 1)];
        num_incor_trials = min([num_incor_trials, size(features(temp_sess_idx).LR, 3), size(features(temp_sess_idx).RL, 3)]);
    end
end
clear results vals vals_cumsum

%% Proportion of direction-selective units
proportions = zeros(numel(cond_name_list), size(param.epochs, 1), param.num_repeat);

for repeat_idx = 1:param.num_repeat
    rng(repeat_idx, 'twister');
    for cond_idx = 1:numel(cond_name_list)
        unit_list = randi(size(selectivity(cond_idx).is_selective, 1), param.num_sampled_units, 1); % with replacement
        proportions(cond_idx, :, repeat_idx) = mean(selectivity(cond_idx).is_selective(unit_list, :), 1);
    end
end

%% Plot proportion of direction-selective units (all epochs)
num_epochs = size(param.epochs, 1);
proportion_mean = mean(proportions, 3);
proportion_sem = std(proportions, 0, 3)/sqrt(param.num_repeat);

% Before vs Atropine: two-sided Mann-Whitney U, Bonferroni across epochs
proportion_u = zeros(num_epochs, 1);
proportion_p_raw = zeros(num_epochs, 1);
for epoch_idx = 1:num_epochs
    x = squeeze(proportions(1, epoch_idx, :));
    y = squeeze(proportions(2, epoch_idx, :));
    ranks = tiedrank([x; y]);
    proportion_u(epoch_idx) = sum(ranks(1:numel(x))) - numel(x)*(numel(x)+1)/2;
    if all([x; y] == x(1))
        proportion_p_raw(epoch_idx) = 1;
    else
        proportion_p_raw(epoch_idx) = ranksum(x, y, 'tail', 'both');
    end
end
proportion_p_bonferroni = min(proportion_p_raw*num_epochs, 1);
proportion_significance = repmat("ns", num_epochs, 1);
proportion_significance(proportion_p_bonferroni < 0.05) = "*";
proportion_significance(proportion_p_bonferroni < 0.01) = "**";
proportion_significance(proportion_p_bonferroni < 0.001) = "***";

stats_proportions = table(param.epoch_names, proportion_u, proportion_p_raw, proportion_p_bonferroni, proportion_significance, ...
    'VariableNames', {'Epoch', 'U', 'P_raw', 'P_Bonferroni', 'Significance'});
disp(stats_proportions);

selected_colors = [0.85, 0.33, 0.10; 0.35, 0.35, 0.35]; % Before, Atropine

% fig = figure('Position', [80, 80, 1250, 720]);
% axes('Position', [0.10, 0.15, 0.86, 0.66]);
% hold on
% h = bar(proportion_mean', 'grouped');
% bar_x = zeros(numel(cond_name_list), num_epochs);
% for cond_idx = 1:numel(cond_name_list)
%     h(cond_idx).FaceColor = selected_colors(cond_idx, :);
%     bar_x(cond_idx, :) = h(cond_idx).XEndPoints;
%     errorbar(bar_x(cond_idx, :), proportion_mean(cond_idx, :), proportion_sem(cond_idx, :), 'k', 'LineStyle', 'none', 'LineWidth', 1.3, 'HandleVisibility', 'off');
% end
% bracket_y = zeros(1, num_epochs);
% for epoch_idx = 1:num_epochs
%     bracket_y(epoch_idx) = max([0; proportion_mean(:, epoch_idx) + proportion_sem(:, epoch_idx)]) + 0.04;
%     temp_x = bar_x(:, epoch_idx);
%     temp_y = bracket_y(epoch_idx);
%     plot([temp_x(1), temp_x(1), temp_x(2), temp_x(2)], [temp_y-0.015, temp_y, temp_y, temp_y-0.015], 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
%     text(mean(temp_x), temp_y+0.006, proportion_significance(epoch_idx), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 13);
% end
% xticks(1:num_epochs);
% xticklabels(param.epoch_names);
% ylim([0, max(1, max(bracket_y)+0.08)]);
% yticks(0:0.2:1);
% ylabel('Proportion of direction-selective units (mean +/- SEM)');
% title({sprintf('%d units per group, %d repeats (with replacement)', param.num_sampled_units, param.num_repeat), ...
%     sprintf('Unit selection: RR vs LL p < %.2f; bars: Mann-Whitney + Bonferroni (%d epochs)', param.sel_alpha, num_epochs)});
% legend(h, cond_name_list, 'Location', 'northoutside', 'Orientation', 'horizontal');
% set(gca, 'FontSize', 13, 'Box', 'off');
% annotation('textbox', [0.10, 0.025, 0.86, 0.04], 'String', '* p < 0.05, ** p < 0.01, *** p < 0.001; ns: not significant (adjusted p)', ...
%     'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 11);
% print(fig, strcat(param.save_path, parts_path, 'EpochDirectionSelectiveProportions'), '-dsvg');

%% Plot proportion of direction-selective units (delay epochs)
% Same p-values as above, Bonferroni across the two delay epochs only
delay_p_bonferroni = min(proportion_p_raw(param.delay_epoch_idx)*numel(param.delay_epoch_idx), 1);
delay_significance = repmat("ns", numel(param.delay_epoch_idx), 1);
delay_significance(delay_p_bonferroni < 0.05) = "*";
delay_significance(delay_p_bonferroni < 0.01) = "**";
delay_significance(delay_p_bonferroni < 0.001) = "***";

stats_proportions_delay = table(param.epoch_names(param.delay_epoch_idx), proportion_u(param.delay_epoch_idx), proportion_p_raw(param.delay_epoch_idx), ...
    delay_p_bonferroni, delay_significance, 'VariableNames', {'Epoch', 'U', 'P_raw', 'P_Bonferroni', 'Significance'});
disp(stats_proportions_delay);

fig = figure('Position', [80, 80, 850, 720]);
axes('Position', [0.10, 0.15, 0.86, 0.66]);
hold on
bar_x = (1:2) + [-1; 1]/7; % same positions as a grouped bar plot
bar_width = 0.8/3.5;
h = gobjects(numel(cond_name_list), 2);
for cond_idx = 1:numel(cond_name_list)
    for plot_idx = 1:2
        temp_mean = proportion_mean(cond_idx, param.delay_epoch_idx(plot_idx));
        h(cond_idx, plot_idx) = patch(bar_x(cond_idx, plot_idx) + [-1, 1, 1, -1]*bar_width/2, [0, 0, temp_mean, temp_mean], selected_colors(cond_idx, :), ...
            'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', 0.5);
    end
    errorbar(bar_x(cond_idx, :), proportion_mean(cond_idx, param.delay_epoch_idx), proportion_sem(cond_idx, param.delay_epoch_idx), 'k', ...
        'LineStyle', 'none', 'LineWidth', 1.3, 'HandleVisibility', 'off');
end
for plot_idx = 1:2
    epoch_idx = param.delay_epoch_idx(plot_idx);
    temp_x = bar_x(:, plot_idx);
    temp_y = max([0; proportion_mean(:, epoch_idx) + proportion_sem(:, epoch_idx)]) + 0.04;
    plot([temp_x(1), temp_x(1), temp_x(2), temp_x(2)], [temp_y-0.015, temp_y, temp_y, temp_y-0.015], 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    text(mean(temp_x), temp_y+0.006, delay_significance(plot_idx), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 13);
end
xlim([0.5, 2.5]);
xticks(1:2);
xticklabels(param.epoch_names(param.delay_epoch_idx));
ylim([0, 0.5]);
yticks(0:0.1:0.5);
ylabel('Proportion of direction-selective units (mean +/- SEM)');
title({sprintf('%d units per group, %d repeats (with replacement)', param.num_sampled_units, param.num_repeat), ...
    sprintf('Unit selection: RR vs LL p < %.2f; bars: Mann-Whitney + Bonferroni (%d epochs)', param.sel_alpha, numel(param.delay_epoch_idx))});
legend(h(:, 1), cond_name_list, 'Location', 'northoutside', 'Orientation', 'horizontal');
set(gca, 'FontSize', 13, 'Box', 'off');
annotation('textbox', [0.10, 0.025, 0.86, 0.04], 'String', '* p < 0.05, ** p < 0.01, *** p < 0.001; ns: not significant (adjusted p)', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 11);
print(fig, strcat(param.save_path, parts_path, 'DelayDirectionSelectiveProportions'), '-dsvg');

%% Plot reduction of direction-selective units by atropine (Delay Early vs Delay Late)
% (P_Before - P_Atropine)/(P_Before - floor) per repeat
% floor = sel_alpha removes the chance level of the selection (fraction of above-chance selective units lost)
floor_list = [0, param.sel_alpha];
fname_list = ["DelaySelectivityReduction"; "DelaySelectivityReduction_ChanceCorrected"];
p_before = reshape(proportions(1, param.delay_epoch_idx, :), 2, []); % epochs x repeats
p_atropine = reshape(proportions(2, param.delay_epoch_idx, :), 2, []);
stats_reduction = struct();

for floor_idx = 1 % :numel(floor_list)
    reduction = (p_before - p_atropine)./(p_before - floor_list(floor_idx));
    reduction(~(p_before > floor_list(floor_idx))) = NaN; % undefined unless Before exceeds the floor
    num_valid = sum(isfinite(reduction), 2);
    reduction_mean = mean(reduction, 2, 'omitnan');
    reduction_sem = std(reduction, 0, 2, 'omitnan')./sqrt(num_valid);

    % Delay Early vs Delay Late: two-sided Mann-Whitney U (single comparison)
    x = reduction(1, isfinite(reduction(1, :)));
    y = reduction(2, isfinite(reduction(2, :)));
    ranks = tiedrank([x(:); y(:)]);
    reduction_u = sum(ranks(1:numel(x))) - numel(x)*(numel(x)+1)/2;
    if all([x(:); y(:)] == x(1))
        reduction_p = 1;
    else
        reduction_p = ranksum(x, y, 'tail', 'both');
    end
    reduction_significance = "ns";
    reduction_significance(reduction_p < 0.05) = "*";
    reduction_significance(reduction_p < 0.01) = "**";
    reduction_significance(reduction_p < 0.001) = "***";

    stats_reduction(floor_idx).floor = floor_list(floor_idx);
    stats_reduction(floor_idx).values = reduction;
    stats_reduction(floor_idx).num_valid = num_valid;
    stats_reduction(floor_idx).mean = reduction_mean;
    stats_reduction(floor_idx).sem = reduction_sem;
    stats_reduction(floor_idx).u = reduction_u;
    stats_reduction(floor_idx).p_raw = reduction_p;
    stats_reduction(floor_idx).significance = reduction_significance;
    fprintf('Reduction (floor %.2f): Delay Early %.4f +/- %.4f, Delay Late %.4f +/- %.4f, U = %g, p = %.4g (%s)\n', ...
        floor_list(floor_idx), reduction_mean(1), reduction_sem(1), reduction_mean(2), reduction_sem(2), reduction_u, reduction_p, reduction_significance);

    fig = figure('Position', [140, 80, 700, 720]);
    axes('Position', [0.17, 0.15, 0.78, 0.64]);
    hold on
    bar_colors = [0.6, 0.6, 0.6; 0.3, 0.3, 0.3];
    bar_width = 0.6;
    for plot_idx = 1:2
        patch(plot_idx + [-1, 1, 1, -1]*bar_width/2, [0, 0, reduction_mean(plot_idx), reduction_mean(plot_idx)], bar_colors(plot_idx, :), ...
            'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', 0.5);
    end
    errorbar(1:2, reduction_mean, reduction_sem, 'k', 'LineStyle', 'none', 'LineWidth', 1.3);
    bracket_y = max([0; reduction_mean + reduction_sem]) + 0.06;
    plot([1, 1, 2, 2], [bracket_y-0.02, bracket_y, bracket_y, bracket_y-0.02], 'k-', 'LineWidth', 1.2);
    text(1.5, bracket_y+0.008, reduction_significance, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 13);
    y_low = min([0; reduction_mean - reduction_sem]);
    if y_low < 0
        yline(0, '-k');
        y_low = y_low - 0.05;
    end
    xlim([0.4, 2.6]);
    ylim([y_low, max(1, bracket_y+0.1)]);
    yticks(-1:0.2:1);
    xticks(1:2);
    xticklabels(param.epoch_names(param.delay_epoch_idx));
    if floor_list(floor_idx) > 0
        ylabel({sprintf('(P_{Before} - P_{Atropine}) / (P_{Before} - %.2f)', floor_list(floor_idx)), '(mean +/- SEM)'});
        title({'Atropine reduction of above-chance selective units', ...
            sprintf('Floor %.2f = false-positive rate of p < %.2f selection', floor_list(floor_idx), param.sel_alpha), ...
            sprintf('%d repeats; two-sided Mann-Whitney U', param.num_repeat)});
    else
        ylabel({'(P_{Before} - P_{Atropine}) / P_{Before}', '(mean +/- SEM)'});
        title({'Atropine reduction of direction-selective units', sprintf('%d repeats; two-sided Mann-Whitney U', param.num_repeat)});
    end
    set(gca, 'FontSize', 13, 'Box', 'off');
    annotation('textbox', [0.05, 0.02, 0.90, 0.06], 'String', {'* p < 0.05, ** p < 0.01, *** p < 0.001; ns: not significant (single test)', ...
        sprintf('Valid repeats: Delay Early %d, Delay Late %d', num_valid)}, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 11);
    print(fig, strcat(param.save_path, parts_path, fname_list(floor_idx)), '-dsvg');
end

