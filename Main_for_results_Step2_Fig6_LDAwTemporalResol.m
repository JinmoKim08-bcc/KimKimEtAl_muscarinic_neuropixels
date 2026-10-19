clc
clear
close all

%% Hyperparameters to adjust
param = struct();

% add path NPY library for reading Phy
npy_path = "utils\npy-matlab\npy-matlab\";
addpath(npy_path)
addpath('utils\hatchfill2_r8')

save_path = 'C:\Users\sirom\Documents\KimKimEtAl_muscarnic_neuropixels\';
parts_path = '\fig6\';

param.version = "v16fd";

fig_path = strcat(save_path, param.version, parts_path);
if ~isfolder(fig_path)
    mkdir(fig_path);
end

%% Hyperparameters that most do not need to be adjusted
param.save_path = save_path;
param.fs = 30000; % [Hz]; sampling rate
param.bin_time = 0.001; % [sec]; temporal binning size to calculate firing rate
param.smoothing_pts = 100; % [idx]; smoothing data points
param.smoothing_window = ones(1, param.smoothing_pts) / param.smoothing_pts;
param.ITI_time = 2; % [sec]; Inter-trial interval duration
param.sample_time = 0.1; % [sec]; timing given a sample
param.delay_time = 2; % [sec]
param.gocue_time = param.sample_time + param.delay_time; % [sec]
param.response_time = 0.5; % [sec]
param.t_max = param.gocue_time + param.response_time; % [sec]; maximal end time to analyze
param.t_max_truncated = param.t_max; % [sec]

% parameters for time-resolved LDA
param.lda_win_time_ms = 200; % [ms]; causal window length
param.lda_step_time_ms = 50; % [ms]
param.start_offset_ms = 1000; % [ms]; first window ends at -0.8 s
param.num_repeat = 100;
param.num_chosen_units = 100;
param.num_trials_per_class = 50; % correct trials (RR / LL) per unit
param.num_folds = 5;
param.min_train_trials = 2;
param.min_test_trials = 1;
param.rng_repeat = 1e4; % seed offset for unit/trial sampling

% parameters for response-selective units
param.sel_alpha = 0.05;

% parameters for bar graphs
param.bar_times = [0.8, 1.5]; % [sec]
param.bar_ylim = [0.4, 1.0];
param.fig_ylim = [0.3, 1.0];

% parameters for population-vector correlation
param.pv_num_repeat = 100;
param.pv_num_units = 100;
param.pv_epoch_delay_early = [0.5, 1.0]; % [sec]
param.pv_epoch_delay_late = [1.6, 2.1]; % [sec]
param.pv_epoch_response = [2.1, 2.6]; % [sec]
param.pv_ylim = [-0.2, 1.0];
param.pv_ylim_primary = [0, 1.0];

cond_name_list = ["Before"; "Atropine"; "Saline"];
plot_cond_list = [1, 2]; % conditions shown in the figures (Before, Atropine)
num_sess_list = [5, 4, 3];
trial_types = ["RR"; "LL"; "RL"; "LR"];
num_cond = numel(cond_name_list);

%% Causal time windows
ITI_dur_bin = round(param.ITI_time/param.bin_time);
num_pts = round((param.ITI_time+param.t_max)/param.bin_time);
time_bin = -param.ITI_time:param.bin_time:(param.t_max-param.bin_time);

win_pts = round(param.lda_win_time_ms/(1000*param.bin_time));
step_pts = round(param.lda_step_time_ms/(1000*param.bin_time));
start_pts = round(param.start_offset_ms/(1000*param.bin_time));
ub_list = (win_pts+start_pts):step_pts:num_pts;
lb_list = ub_list - win_pts + 1;
time_bin_win = time_bin(ub_list);
num_win = numel(ub_list);

%% Read calculated firing vectors
load(strcat(param.save_path, 'results_', param.version, '.mat'), 'results');

unit_info = struct();
for cond_idx = 1:num_cond
    num_units_list = zeros(num_sess_list(cond_idx), 1);
    num_trials_list = zeros(num_sess_list(cond_idx), numel(trial_types));
    for sess_idx = 1:num_sess_list(cond_idx)
        num_units_list(sess_idx) = size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1);
        for trial_idx = 1:numel(trial_types)
            num_trials_list(sess_idx, trial_idx) = size(results(cond_idx).sess(sess_idx).(trial_types(trial_idx)).allTrials.vals, 3);
        end
    end
    num_units_cumsum2 = cumsum(num_units_list);
    num_units_cumsum1 = [1; num_units_cumsum2(1:end-1)+1];

    unit_info(cond_idx).num_units_list = num_units_list;
    unit_info(cond_idx).num_trials_list = num_trials_list;
    unit_info(cond_idx).num_units_cumsum = [num_units_cumsum1, num_units_cumsum2];
    unit_info(cond_idx).num_units = num_units_cumsum2(end);
end

% incorrect trials: use the trial count of the sparsest session (RL / LR)
min_incorrect = inf;
for cond_idx = 1:num_cond
    min_incorrect = min(min_incorrect, min(min(unit_info(cond_idx).num_trials_list(:, 3:4))));
end

% class R is the right choice in both pairs (LR = stimulus L, answered R)
pair_info = struct();
pair_info(1).name = 'Correct: RR vs LL';
pair_info(1).short = 'correct';
pair_info(1).type_R = 'RR';
pair_info(1).type_L = 'LL';
pair_info(1).num_trials = param.num_trials_per_class;
pair_info(1).num_folds = param.num_folds;

pair_info(2).name = 'Incorrect: LR vs RL';
pair_info(2).short = 'incorrect';
pair_info(2).type_R = 'LR';
pair_info(2).type_L = 'RL';
pair_info(2).num_trials = min_incorrect;
pair_info(2).num_folds = min(param.num_folds, min_incorrect);

num_pair = numel(pair_info);

%% Response-selective units (statistics): RR vs LL ranksum in [2.1, 2.6] s
resp_lb_idx = round(param.gocue_time/param.bin_time)+1;
resp_ub_idx = round(param.t_max/param.bin_time);

sel_stat = struct();
for cond_idx = 1:num_cond
    sel_idx_pooled = [];
    for sess_idx = 1:num_sess_list(cond_idx)
        epoch_fr = cell(2, 1);
        for trial_idx = 1:2
            vals = results(cond_idx).sess(sess_idx).(trial_types(trial_idx)).allTrials.vals;

            baseline = vals(:, 1:ITI_dur_bin, :);
            baseline(isnan(baseline)) = 0;
            baseline = mean(baseline, 2);

            temp_fr = vals(:, ITI_dur_bin+1:num_pts, :);
            temp_fr(isnan(temp_fr)) = 0;
            temp_fr = temp_fr - baseline;
            temp_fr = movmean(temp_fr, [param.smoothing_pts-1, 0], 2, 'Endpoints', 'shrink'); % causal smoothing

            epoch_fr{trial_idx} = reshape(mean(temp_fr(:, resp_lb_idx:resp_ub_idx, :), 2), size(temp_fr, 1), []); % unit x trial
        end

        num_units = size(epoch_fr{1}, 1);
        p_val = ones(num_units, 1);
        for unit_idx = 1:num_units
            x = epoch_fr{1}(unit_idx, :);
            y = epoch_fr{2}(unit_idx, :);
            if all([x, y] == x(1))
                continue
            end
            p_val(unit_idx) = ranksum(x, y);
        end

        sel_stat(cond_idx).sess(sess_idx).p_val = p_val;
        sel_stat(cond_idx).sess(sess_idx).is_sel = p_val < param.sel_alpha;
        sel_idx_pooled = [sel_idx_pooled; unit_info(cond_idx).num_units_cumsum(sess_idx, 1)-1+find(p_val < param.sel_alpha)];

        fprintf('%s, Sess %d: %d / %d response-selective units\n', cond_name_list(cond_idx), sess_idx, sum(p_val < param.sel_alpha), num_units);
    end
    unit_info(cond_idx).sel_stat_idx = sel_idx_pooled;
end

%% Populations: full pool vs response-selective units removed
pop_info = struct();
pop_info(1).name = 'Full';
pop_info(1).label = 'Full (intact)';
pop_info(2).name = 'dropSel_Response';
pop_info(2).label = 'Response-selective removed [statistics, 2.1-2.6 s]';
num_pop = numel(pop_info);

pool_sizes = zeros(num_pop, num_cond);
for cond_idx = 1:num_cond
    drop_idx_list = {[], unit_info(cond_idx).sel_stat_idx};
    for pop_idx = 1:num_pop
        drop_idx = unique(drop_idx_list{pop_idx});
        pop_info(pop_idx).drop_idx{cond_idx} = drop_idx;
        pop_info(pop_idx).pool_idx{cond_idx} = setdiff((1:unit_info(cond_idx).num_units)', drop_idx);
        pool_sizes(pop_idx, cond_idx) = numel(pop_info(pop_idx).pool_idx{cond_idx});
    end
end
disp(pool_sizes)

% same number of units for every population and condition
num_units_eff = min(param.num_chosen_units, min(pool_sizes(:)));

%% Causal single-trial features
% win_feat{cond, sess}.(trial type): unit x window x trial, baseline-subtracted per trial
win_feat = cell(num_cond, max(num_sess_list));
for cond_idx = 1:num_cond
    for sess_idx = 1:num_sess_list(cond_idx)
        win_feat{cond_idx, sess_idx} = struct();
        for trial_idx = 1:numel(trial_types)
            selected_trial_type = trial_types(trial_idx);

            vals = single(results(cond_idx).sess(sess_idx).(selected_trial_type).allTrials.vals(:, 1:num_pts, :));
            vals(isnan(vals)) = 0;
            baseline = mean(vals(:, 1:ITI_dur_bin, :), 2);
            vals_cumsum = cumsum(vals, 2);

            temp_feat = zeros(size(vals, 1), num_win, size(vals, 3), 'single');
            for win_idx = 1:num_win
                temp_mean = (vals_cumsum(:, ub_list(win_idx), :) - vals_cumsum(:, lb_list(win_idx)-1, :)) / single(ub_list(win_idx)-lb_list(win_idx)+1);
                temp_feat(:, win_idx, :) = temp_mean - baseline;
            end
            win_feat{cond_idx, sess_idx}.(selected_trial_type) = temp_feat;
        end
    end
    fprintf('%s: features done\n', cond_name_list(cond_idx));
end

%% Population vectors for the correlation analysis
% trial-averaged, smoothed and z-scored to the ITI (same as Fig5)
pv_time_bin = -1:param.bin_time:(param.t_max_truncated-param.bin_time);
pv_keep_from = round(1/param.bin_time)+1; % drop the first second so that the time axis starts at -1 s

pv = struct();
for cond_idx = 1:num_cond
    CR_all = [];
    CL_all = [];
    ER_all = [];
    EL_all = [];
    for sess_idx = 1:num_sess_list(cond_idx)
        temp_CR = results(cond_idx).sess(sess_idx).RR.allTrials.vals(:, 1:num_pts, :);
        temp_CL = results(cond_idx).sess(sess_idx).LL.allTrials.vals(:, 1:num_pts, :);
        temp_ER = results(cond_idx).sess(sess_idx).LR.allTrials.vals(:, 1:num_pts, :);
        temp_EL = results(cond_idx).sess(sess_idx).RL.allTrials.vals(:, 1:num_pts, :);
        temp_results = Calculate_all_even_odd_trials_vector(temp_CR, temp_CL, temp_ER, temp_EL, param);

        CR_all = [CR_all; temp_results.CR_all];
        CL_all = [CL_all; temp_results.CL_all];
        ER_all = [ER_all; temp_results.ER_all];
        EL_all = [EL_all; temp_results.EL_all];
    end
    pv(cond_idx).CR = CR_all(:, pv_keep_from:end);
    pv(cond_idx).CL = CL_all(:, pv_keep_from:end);
    pv(cond_idx).ER = ER_all(:, pv_keep_from:end);
    pv(cond_idx).EL = EL_all(:, pv_keep_from:end);
end

clear results temp_CR temp_CL temp_ER temp_EL temp_results CR_all CL_all ER_all EL_all

%% LDA analysis with time resolution (matched trials, k-fold CV)
lda_results = struct();

for pop_idx = 1:num_pop
    for cond_idx = 1:num_cond
        num_units_cumsum = unit_info(cond_idx).num_units_cumsum;
        pool_idx = pop_info(pop_idx).pool_idx{cond_idx};
        seed_base = param.rng_repeat + cond_idx*1000;

        auroc = nan(num_win, param.num_repeat, num_pair);
        within_var = nan(num_win, param.num_repeat, num_pair);
        selected_units = nan(num_units_eff, param.num_repeat);

        for repeat_idx = 1:param.num_repeat
            fprintf('Pop: [%d/%d], Cond: [%d/%d], Repeat: [%d/%d]\n', pop_idx, num_pop, cond_idx, num_cond, repeat_idx, param.num_repeat);

            % sample units
            rng(seed_base+repeat_idx, 'twister');
            chosen_units = pool_idx(randperm(numel(pool_idx), num_units_eff))';
            selected_units(:, repeat_idx) = chosen_units(:);

            % session and row of each sampled unit
            unit_sess = zeros(num_units_eff, 1);
            unit_row = zeros(num_units_eff, 1);
            for unit_idx = 1:num_units_eff
                unit_sess(unit_idx) = find(chosen_units(unit_idx) >= num_units_cumsum(:, 1) & chosen_units(unit_idx) <= num_units_cumsum(:, 2), 1, 'first');
                unit_row(unit_idx) = chosen_units(unit_idx) - num_units_cumsum(unit_sess(unit_idx), 1) + 1;
            end

            for pair_idx = 1:num_pair
                rng(seed_base+repeat_idx*10+pair_idx, 'twister');
                type_R = pair_info(pair_idx).type_R;
                type_L = pair_info(pair_idx).type_L;
                num_trials = pair_info(pair_idx).num_trials;
                num_folds = pair_info(pair_idx).num_folds;

                % fold sizes differ by at most one trial
                fold_size = floor(num_trials/num_folds)*ones(1, num_folds);
                fold_size(1:mod(num_trials, num_folds)) = fold_size(1:mod(num_trials, num_folds)) + 1;
                fold_template = repelem(1:num_folds, fold_size);

                % same number of R and L trials for every unit, split into folds
                trial_R = cell(num_units_eff, 1);
                trial_L = cell(num_units_eff, 1);
                fold_R = cell(num_units_eff, 1);
                fold_L = cell(num_units_eff, 1);
                for unit_idx = 1:num_units_eff
                    trial_R{unit_idx} = randperm(size(win_feat{cond_idx, unit_sess(unit_idx)}.(type_R), 3), num_trials);
                    trial_L{unit_idx} = randperm(size(win_feat{cond_idx, unit_sess(unit_idx)}.(type_L), 3), num_trials);

                    temp_fold = zeros(1, num_trials);
                    temp_fold(randperm(num_trials)) = fold_template;
                    fold_R{unit_idx} = temp_fold;

                    temp_fold = zeros(1, num_trials);
                    temp_fold(randperm(num_trials)) = fold_template;
                    fold_L{unit_idx} = temp_fold;
                end

                fold_auroc = nan(num_win, param.num_folds);
                fold_var = nan(num_win, param.num_folds);
                fold_weight = zeros(1, param.num_folds);

                for fold_idx = 1:num_folds
                    num_test = fold_size(fold_idx);
                    num_train = num_trials - num_test;

                    % pseudo population: trials are drawn independently for each unit
                    X_train_R = zeros(num_train, num_units_eff, num_win, 'single');
                    X_train_L = zeros(num_train, num_units_eff, num_win, 'single');
                    X_test_R = zeros(num_test, num_units_eff, num_win, 'single');
                    X_test_L = zeros(num_test, num_units_eff, num_win, 'single');

                    for unit_idx = 1:num_units_eff
                        temp_trials = trial_R{unit_idx}(fold_R{unit_idx} ~= fold_idx);
                        temp_trials = temp_trials(randperm(numel(temp_trials), num_train));
                        X_train_R(:, unit_idx, :) = permute(win_feat{cond_idx, unit_sess(unit_idx)}.(type_R)(unit_row(unit_idx), :, temp_trials), [3, 1, 2]);
                    end
                    for unit_idx = 1:num_units_eff
                        temp_trials = trial_L{unit_idx}(fold_L{unit_idx} ~= fold_idx);
                        temp_trials = temp_trials(randperm(numel(temp_trials), num_train));
                        X_train_L(:, unit_idx, :) = permute(win_feat{cond_idx, unit_sess(unit_idx)}.(type_L)(unit_row(unit_idx), :, temp_trials), [3, 1, 2]);
                    end
                    for unit_idx = 1:num_units_eff
                        temp_trials = trial_R{unit_idx}(fold_R{unit_idx} == fold_idx);
                        temp_trials = temp_trials(randperm(numel(temp_trials), num_test));
                        X_test_R(:, unit_idx, :) = permute(win_feat{cond_idx, unit_sess(unit_idx)}.(type_R)(unit_row(unit_idx), :, temp_trials), [3, 1, 2]);
                    end
                    for unit_idx = 1:num_units_eff
                        temp_trials = trial_L{unit_idx}(fold_L{unit_idx} == fold_idx);
                        temp_trials = temp_trials(randperm(numel(temp_trials), num_test));
                        X_test_L(:, unit_idx, :) = permute(win_feat{cond_idx, unit_sess(unit_idx)}.(type_L)(unit_row(unit_idx), :, temp_trials), [3, 1, 2]);
                    end

                    if num_train < param.min_train_trials || num_test < param.min_test_trials
                        continue
                    end
                    fold_weight(fold_idx) = num_test^2;

                    for win_idx = 1:num_win
                        X_train = double([X_train_R(:, :, win_idx); X_train_L(:, :, win_idx)]);
                        X_test = double([X_test_R(:, :, win_idx); X_test_L(:, :, win_idx)]);

                        % pseudoLinear LDA (same projection as fitcdiscr 'pseudoLinear')
                        mu_train = mean(X_train, 1);
                        X_train = X_train - mu_train;
                        X_test = X_test - mu_train;

                        X_R = X_train(1:num_train, :);
                        X_L = X_train(num_train+1:end, :);
                        m_R = mean(X_R, 1);
                        m_L = mean(X_L, 1);
                        C_R = X_R - m_R;
                        C_L = X_L - m_L;
                        Sw = (C_R'*C_R + C_L'*C_L) / (2*num_train-2);
                        fold_var(win_idx, fold_idx) = trace(Sw)/size(Sw, 1);

                        w = pinv(Sw)*(m_R-m_L)';
                        if all(isfinite(w)) && norm(w) > 0
                            score_rank = tiedrank(X_test*w);
                            fold_auroc(win_idx, fold_idx) = (sum(score_rank(1:num_test)) - num_test*(num_test+1)/2) / (num_test*num_test);
                        end
                    end
                end

                % weighted mean over folds (weight = num_test^2)
                weight_mat = repmat(fold_weight, num_win, 1) .* isfinite(fold_auroc);
                temp_auroc = fold_auroc;
                temp_auroc(weight_mat == 0) = 0;
                auroc(:, repeat_idx, pair_idx) = sum(temp_auroc.*weight_mat, 2) ./ sum(weight_mat, 2);

                weight_mat = repmat(fold_weight, num_win, 1) .* isfinite(fold_var);
                temp_var = fold_var;
                temp_var(weight_mat == 0) = 0;
                within_var(:, repeat_idx, pair_idx) = sum(temp_var.*weight_mat, 2) ./ sum(weight_mat, 2);
            end
        end

        lda_results(pop_idx).cond(cond_idx).auroc = auroc;
        lda_results(pop_idx).cond(cond_idx).within_var = within_var;
        lda_results(pop_idx).cond(cond_idx).selected_units = selected_units;
    end
end

save(strcat(param.save_path, 'lda_results_', param.version, '-Fig6CandF-matched_trials_Unit', num2str(num_units_eff), '.mat'), 'lda_results', 'pop_info', 'pair_info', 'time_bin_win', 'num_units_eff', '-v7.3')
load(strcat(param.save_path, 'lda_results_', param.version, '-Fig6CandF-matched_trials_Unit', num2str(num_units_eff), '.mat'))

%% Print auROC at a few time points
report_times = [0.8, 1.5, 2.0, 2.3, 2.5];
for pop_idx = 1:num_pop
    fprintf('\nPopulation: %s\n', pop_info(pop_idx).name);
    for cond_idx = 1:num_cond
        fprintf('  %s\n', cond_name_list(cond_idx));
        for time_idx = 1:numel(report_times)
            [~, win_idx] = min(abs(time_bin_win - report_times(time_idx)));
            fprintf('    %.3f s: correct %.3f, incorrect %.3f\n', time_bin_win(win_idx), ...
                mean(lda_results(pop_idx).cond(cond_idx).auroc(win_idx, :, 1), 'omitnan'), ...
                mean(lda_results(pop_idx).cond(cond_idx).auroc(win_idx, :, 2), 'omitnan'));
        end
    end
end

%% Plot time-resolved auROC (all populations)
cond_line = [0.00 0.65 0.10; 1.00 0.55 0.00; 0.49 0.18 0.56];
cond_fill = [0.50 1.00 0.50; 1.00 0.75 0.45; 0.72 0.55 0.78];

for pair_idx = 1:num_pair
    fig = figure('Position', [40 40 min(1860, 540*num_pop) 800]);
    for pop_idx = 1:num_pop
        %%% auROC
        subplot(2, num_pop, pop_idx)
        hold on
        h_line = gobjects(numel(plot_cond_list), 1);
        for plot_idx = 1:numel(plot_cond_list)
            cond_idx = plot_cond_list(plot_idx);
            temp_auroc = lda_results(pop_idx).cond(cond_idx).auroc(:, :, pair_idx);
            meanAUC = mean(temp_auroc, 2, 'omitnan');
            semAUC = std(temp_auroc, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(temp_auroc), 2));

            h = fill([time_bin_win'; flipud(time_bin_win')], [meanAUC+semAUC; flipud(meanAUC-semAUC)], cond_fill(cond_idx, :), 'LineStyle', 'none', 'HandleVisibility', 'off');
            set(h, 'FaceAlpha', 0.25);
            h_line(plot_idx) = plot(time_bin_win, meanAUC, 'Color', cond_line(cond_idx, :), 'LineWidth', 2.8);

            % Full population as a dashed reference
            if pop_idx ~= 1
                plot(time_bin_win, mean(lda_results(1).cond(cond_idx).auroc(:, :, pair_idx), 2, 'omitnan'), '--', 'Color', cond_line(cond_idx, :), 'LineWidth', 1.2, 'HandleVisibility', 'off');
            end
        end
        yline(0.5, 'k--', 'LineWidth', 0.8);
        xline(0, ':k');
        xline(param.sample_time, ':k');
        xline(param.gocue_time, ':k', 'LineWidth', 1.3);
        for time_idx = 1:numel(param.bar_times)
            xline(param.bar_times(time_idx), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        xlim([time_bin_win(1), time_bin_win(end)]);
        ylim(param.fig_ylim);
        grid on
        xlabel('Causal window end time [s]');
        ylabel('auROC');
        title(pop_info(pop_idx).label, 'Interpreter', 'none');
        legend(h_line, cellstr(cond_name_list(plot_cond_list)), 'Location', 'southeast');
        set(gca, 'FontSize', 11);

        %%% within-class variance of the training set
        subplot(2, num_pop, num_pop+pop_idx)
        hold on
        for cond_idx = plot_cond_list
            temp_var = lda_results(pop_idx).cond(cond_idx).within_var(:, :, pair_idx);
            meanVar = mean(temp_var, 2, 'omitnan');
            semVar = std(temp_var, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(temp_var), 2));

            h = fill([time_bin_win'; flipud(time_bin_win')], [meanVar+semVar; flipud(meanVar-semVar)], cond_fill(cond_idx, :), 'LineStyle', 'none', 'HandleVisibility', 'off');
            set(h, 'FaceAlpha', 0.25);
            plot(time_bin_win, meanVar, 'Color', cond_line(cond_idx, :), 'LineWidth', 2.8);
        end
        xline(0, ':k');
        xline(param.sample_time, ':k');
        xline(param.gocue_time, ':k', 'LineWidth', 1.3);
        xlim([time_bin_win(1), time_bin_win(end)]);
        grid on
        xlabel('Causal window end time [s]');
        ylabel('trace(S_W)/p');
        title('Training covariance scale');
        set(gca, 'FontSize', 11);
    end
    sgtitle(sprintf('%s | %d-fold CV | %d trials/class | %d units | R=%d | dashed = Full', pair_info(pair_idx).name, pair_info(pair_idx).num_folds, pair_info(pair_idx).num_trials, num_units_eff, param.num_repeat), 'Interpreter', 'none');
    print(fig, strcat(fig_path, 'Fig6CandF-', pair_info(pair_idx).short, '-timeseries'), '-dsvg');
end

%% Bar graph of auROC at 0.8 s vs 1.5 s
% repeats are resamples of the same recordings, not independent animals
bar_idx = zeros(1, numel(param.bar_times));
for time_idx = 1:numel(param.bar_times)
    [~, bar_idx(time_idx)] = min(abs(time_bin_win - param.bar_times(time_idx)));
end
bar_time_actual = time_bin_win(bar_idx);

cond_pairs = [1 3; 1 2; 3 2]; % Before-Saline, Before-Atropine, Saline-Atropine
num_plot_cond = numel(plot_cond_list);
bar_labels = cell(1, 2*num_plot_cond);
for plot_idx = 1:num_plot_cond
    bar_labels{2*plot_idx-1} = sprintf('%s-%.1fs', cond_name_list(plot_cond_list(plot_idx)), bar_time_actual(1));
    bar_labels{2*plot_idx} = sprintf('%s-%.1fs', cond_name_list(plot_cond_list(plot_idx)), bar_time_actual(2));
end

stat_bar = cell(num_pop, num_pair);
for pair_idx = 1:num_pair
    fig = figure('Position', [60 60 min(1800, 560*num_pop) 520]);
    fprintf('\n%s: auROC at %.3f s vs %.3f s\n', pair_info(pair_idx).name, bar_time_actual(1), bar_time_actual(2));

    for pop_idx = 1:num_pop
        X = nan(param.num_repeat, num_cond, 2);
        for cond_idx = 1:num_cond
            X(:, cond_idx, 1) = lda_results(pop_idx).cond(cond_idx).auroc(bar_idx(1), :, pair_idx)';
            X(:, cond_idx, 2) = lda_results(pop_idx).cond(cond_idx).auroc(bar_idx(2), :, pair_idx)';
        end
        mean_auroc = squeeze(mean(X, 1, 'omitnan')); % cond x time
        sem_auroc = squeeze(std(X, 0, 1, 'omitnan')) ./ sqrt(squeeze(sum(isfinite(X), 1)));

        subplot(1, num_pop, pop_idx)
        hold on
        Plot_bar_with_hatch(gca, mean_auroc(plot_cond_list, :), sem_auroc(plot_cond_list, :), cond_line(plot_cond_list, :));

        %%% statistics
        p_within = nan(num_cond, 1);
        for cond_idx = 1:num_cond
            p_within(cond_idx) = signrank(X(:, cond_idx, 1), X(:, cond_idx, 2));
        end
        p_kw = nan(1, 2);
        p_between = nan(size(cond_pairs, 1), 2);
        for time_idx = 1:2
            p_kw(time_idx) = kruskalwallis(X(:, :, time_idx), [], 'off');
            for comp_idx = 1:size(cond_pairs, 1)
                p_between(comp_idx, time_idx) = ranksum(X(:, cond_pairs(comp_idx, 1), time_idx), X(:, cond_pairs(comp_idx, 2), time_idx));
            end
        end
        p_between_bonf = min(1, p_between*size(cond_pairs, 1));

        temp_stat = struct();
        temp_stat.mean = mean_auroc;
        temp_stat.sem = sem_auroc;
        temp_stat.p_within = p_within;
        temp_stat.p_kw = p_kw;
        temp_stat.p_between = p_between;
        temp_stat.p_between_bonf = p_between_bonf;
        stat_bar{pop_idx, pair_idx} = temp_stat;

        fprintf('  %s\n', pop_info(pop_idx).name);
        for cond_idx = 1:num_cond
            fprintf('    %s: %.3f +- %.3f vs %.3f +- %.3f, signrank p = %.3g %s\n', cond_name_list(cond_idx), ...
                mean_auroc(cond_idx, 1), sem_auroc(cond_idx, 1), mean_auroc(cond_idx, 2), sem_auroc(cond_idx, 2), ...
                p_within(cond_idx), Convert_pvalue_to_asterisk(p_within(cond_idx)));
        end
        for time_idx = 1:2
            fprintf('    %.3f s: Kruskal-Wallis p = %.3g\n', bar_time_actual(time_idx), p_kw(time_idx));
            for comp_idx = 1:size(cond_pairs, 1)
                fprintf('      %s vs %s: ranksum p = %.3g, Bonferroni p = %.3g %s\n', ...
                    cond_name_list(cond_pairs(comp_idx, 1)), cond_name_list(cond_pairs(comp_idx, 2)), ...
                    p_between(comp_idx, time_idx), p_between_bonf(comp_idx, time_idx), ...
                    Convert_pvalue_to_asterisk(p_between_bonf(comp_idx, time_idx)));
            end
        end

        %%% significance brackets (shown conditions only)
        tick = 0.018*diff(param.bar_ylim);
        temp_top = mean_auroc(plot_cond_list, :) + sem_auroc(plot_cond_list, :);
        y_pos = max(temp_top(:)) + 2*tick;
        for plot_idx = 1:num_plot_cond
            x_pos = 2*(plot_idx-1)+[1, 2];
            Plot_sig_bracket(gca, x_pos(1), x_pos(2), y_pos, p_within(plot_cond_list(plot_idx)), tick, 10);
        end
        y_pos = y_pos + 3.2*tick;
        for time_idx = 1:2
            % Before vs Atropine (cond_pairs row 2, Bonferroni across the 3 pairs)
            x1 = 2*(find(plot_cond_list == 1)-1)+time_idx;
            x2 = 2*(find(plot_cond_list == 2)-1)+time_idx;
            Plot_sig_bracket(gca, x1, x2, y_pos, p_between_bonf(2, time_idx), tick, 10);
            y_pos = y_pos + 3.2*tick;
        end

        xticks(1:2*num_plot_cond);
        xticklabels(bar_labels);
        xtickangle(35);
        xlim([0.3, 2*num_plot_cond+0.7]);
        ylim([param.bar_ylim(1), max(param.bar_ylim(2), y_pos+tick)]);
        yticks(param.bar_ylim(1):0.1:param.bar_ylim(2));
        yline(0.5, 'k--', 'LineWidth', 0.8);
        ylabel('auROC');
        title(pop_info(pop_idx).label, 'Interpreter', 'none');
        set(gca, 'FontSize', 11);
    end
    sgtitle(sprintf('%s | auROC at %.2f s (plain) vs %.2f s (hatched) | %d units | R=%d', pair_info(pair_idx).name, bar_time_actual(1), bar_time_actual(2), num_units_eff, param.num_repeat), 'Interpreter', 'none');
    print(fig, strcat(fig_path, 'Fig6CandF-', pair_info(pair_idx).short, '-bar0_8vs1_5'), '-dsvg');
end

%% Population-vector correlation across neurons (correct trials)
% history: delay late vs delay early, upcoming: delay late vs response
[~, early_lb] = min(abs(pv_time_bin - param.pv_epoch_delay_early(1)));
[~, early_ub] = min(abs(pv_time_bin - param.pv_epoch_delay_early(2)));
early_ub = early_ub - 1;
[~, late_lb] = min(abs(pv_time_bin - param.pv_epoch_delay_late(1)));
[~, late_ub] = min(abs(pv_time_bin - param.pv_epoch_delay_late(2)));
late_ub = late_ub - 1;
[~, resp_lb] = min(abs(pv_time_bin - param.pv_epoch_response(1)));
[~, resp_ub] = min(abs(pv_time_bin - param.pv_epoch_response(2)));

pv_num_stable = zeros(1, num_cond);
for cond_idx = 1:num_cond
    pv_num_stable(cond_idx) = size(pv(cond_idx).CR, 1);
end
pv_units_eff = min(param.pv_num_units, min(pv_num_stable));

rho = struct();
for cond_idx = 1:num_cond
    history_RR = nan(param.pv_num_repeat, 1);
    history_LL = nan(param.pv_num_repeat, 1);
    upcoming_RR = nan(param.pv_num_repeat, 1);
    upcoming_LL = nan(param.pv_num_repeat, 1);

    for repeat_idx = 1:param.pv_num_repeat
        rng(repeat_idx)
        pv_units = randperm(pv_num_stable(cond_idx), pv_units_eff);

        late_R = mean(pv(cond_idx).CR(pv_units, late_lb:late_ub), 2);
        early_R = mean(pv(cond_idx).CR(pv_units, early_lb:early_ub), 2);
        resp_R = mean(pv(cond_idx).CR(pv_units, resp_lb:resp_ub), 2);
        late_L = mean(pv(cond_idx).CL(pv_units, late_lb:late_ub), 2);
        early_L = mean(pv(cond_idx).CL(pv_units, early_lb:early_ub), 2);
        resp_L = mean(pv(cond_idx).CL(pv_units, resp_lb:resp_ub), 2);

        % RR and LL are correlated separately
        history_RR(repeat_idx) = corr(late_R, early_R);
        upcoming_RR(repeat_idx) = corr(late_R, resp_R);
        history_LL(repeat_idx) = corr(late_L, early_L);
        upcoming_LL(repeat_idx) = corr(late_L, resp_L);
    end

    rho(cond_idx).history_RR = history_RR;
    rho(cond_idx).history_LL = history_LL;
    rho(cond_idx).upcoming_RR = upcoming_RR;
    rho(cond_idx).upcoming_LL = upcoming_LL;
    % mean of the RR and LL correlations (same as r_cor in Fig4)
    rho(cond_idx).history = (history_RR+history_LL)/2;
    rho(cond_idx).upcoming = (upcoming_RR+upcoming_LL)/2;

    fprintf('%s: history r = %.3f, upcoming r = %.3f\n', cond_name_list(cond_idx), mean(rho(cond_idx).history), mean(rho(cond_idx).upcoming));
end

%%% statistics for RR, LL, and mean of RR and LL
pv_measure_list = {'RR', 'LL', 'mean_RR_LL'};
pv_measure_label = {'RR only', 'LL only', 'mean of RR and LL r'};
num_measure = numel(pv_measure_list);

stat_pv = struct();
for measure_idx = 1:num_measure
    X = nan(param.pv_num_repeat, num_cond, 2);
    for cond_idx = 1:num_cond
        switch pv_measure_list{measure_idx}
            case 'RR'
                X(:, cond_idx, 1) = rho(cond_idx).history_RR;
                X(:, cond_idx, 2) = rho(cond_idx).upcoming_RR;
            case 'LL'
                X(:, cond_idx, 1) = rho(cond_idx).history_LL;
                X(:, cond_idx, 2) = rho(cond_idx).upcoming_LL;
            case 'mean_RR_LL'
                X(:, cond_idx, 1) = rho(cond_idx).history;
                X(:, cond_idx, 2) = rho(cond_idx).upcoming;
        end
    end

    p_within = nan(num_cond, 1);
    for cond_idx = 1:num_cond
        p_within(cond_idx) = signrank(X(:, cond_idx, 1), X(:, cond_idx, 2));
    end
    p_kw = nan(1, 2);
    p_between = nan(size(cond_pairs, 1), 2);
    for time_idx = 1:2
        p_kw(time_idx) = kruskalwallis(X(:, :, time_idx), [], 'off');
        for comp_idx = 1:size(cond_pairs, 1)
            p_between(comp_idx, time_idx) = ranksum(X(:, cond_pairs(comp_idx, 1), time_idx), X(:, cond_pairs(comp_idx, 2), time_idx));
        end
    end

    stat_pv(measure_idx).measure = pv_measure_list{measure_idx};
    stat_pv(measure_idx).mean = squeeze(mean(X, 1, 'omitnan'));
    stat_pv(measure_idx).sem = squeeze(std(X, 0, 1, 'omitnan')) ./ sqrt(param.pv_num_repeat);
    stat_pv(measure_idx).p_within = p_within;
    stat_pv(measure_idx).p_kw = p_kw;
    stat_pv(measure_idx).p_between = p_between;
    stat_pv(measure_idx).p_between_bonf = min(1, p_between*size(cond_pairs, 1));

    fprintf('\n%s\n', pv_measure_list{measure_idx});
    for cond_idx = 1:num_cond
        fprintf('  %s: history %.3f +- %.3f, upcoming %.3f +- %.3f, signrank p = %.3g %s\n', cond_name_list(cond_idx), ...
            stat_pv(measure_idx).mean(cond_idx, 1), stat_pv(measure_idx).sem(cond_idx, 1), ...
            stat_pv(measure_idx).mean(cond_idx, 2), stat_pv(measure_idx).sem(cond_idx, 2), ...
            p_within(cond_idx), Convert_pvalue_to_asterisk(p_within(cond_idx)));
    end
    for time_idx = 1:2
        fprintf('  Kruskal-Wallis p = %.3g\n', p_kw(time_idx));
        for comp_idx = 1:size(cond_pairs, 1)
            fprintf('    %s vs %s: ranksum p = %.3g, Bonferroni p = %.3g %s\n', ...
                cond_name_list(cond_pairs(comp_idx, 1)), cond_name_list(cond_pairs(comp_idx, 2)), ...
                p_between(comp_idx, time_idx), stat_pv(measure_idx).p_between_bonf(comp_idx, time_idx), ...
                Convert_pvalue_to_asterisk(stat_pv(measure_idx).p_between_bonf(comp_idx, time_idx)));
        end
    end
end

pv_labels = cell(1, 2*num_plot_cond);
for plot_idx = 1:num_plot_cond
    pv_labels{2*plot_idx-1} = sprintf('%s-History', cond_name_list(plot_cond_list(plot_idx)));
    pv_labels{2*plot_idx} = sprintf('%s-Upcoming', cond_name_list(plot_cond_list(plot_idx)));
end
pv_title = sprintf('Correct trials | delay late [%.1f,%.1f] vs delay early [%.1f,%.1f] (plain) and response [%.1f,%.1f] (hatched) | %d units x %d repeats', ...
    param.pv_epoch_delay_late, param.pv_epoch_delay_early, param.pv_epoch_response, pv_units_eff, param.pv_num_repeat);

%%% RR, LL, and mean of RR and LL side by side
fig = figure('Position', [60 60 min(1860, 620*num_measure) 560]);
for measure_idx = 1:num_measure
    subplot(1, num_measure, measure_idx)
    hold on
    Plot_bar_with_hatch(gca, stat_pv(measure_idx).mean(plot_cond_list, :), stat_pv(measure_idx).sem(plot_cond_list, :), cond_line(plot_cond_list, :));
    xticks(1:2*num_plot_cond);
    xticklabels(pv_labels);
    xtickangle(35);
    xlim([0.3, 2*num_plot_cond+0.7]);
    ylim(param.pv_ylim);
    yticks(param.pv_ylim(1):0.2:param.pv_ylim(2));
    yline(0, 'k-', 'LineWidth', 0.8);
    ylabel('Correlation coefficient across neurons');
    title({pv_measure_list{measure_idx}, pv_measure_label{measure_idx}}, 'Interpreter', 'none', 'FontSize', 10);
    set(gca, 'FontSize', 11);
end
sgtitle(pv_title, 'Interpreter', 'none');
print(fig, strcat(fig_path, 'Fig6CandF-PopVecCorr-history-vs-upcoming-RR-LL-mean_RR_LL'), '-dsvg');

%%% mean of RR and LL only
fig = figure('Position', [100 100 720 560]);
hold on
Plot_bar_with_hatch(gca, stat_pv(3).mean(plot_cond_list, :), stat_pv(3).sem(plot_cond_list, :), cond_line(plot_cond_list, :));
xticks(1:2*num_plot_cond);
xticklabels(pv_labels);
xtickangle(35);
xlim([0.3, 2*num_plot_cond+0.7]);
ylim(param.pv_ylim_primary);
yticks(param.pv_ylim_primary(1):0.2:param.pv_ylim_primary(2));
yline(0, 'k-', 'LineWidth', 0.8);
ylabel('Correlation coefficient across neurons');
title({pv_measure_list{3}, pv_measure_label{3}}, 'Interpreter', 'none', 'FontSize', 10);
set(gca, 'FontSize', 11);
sgtitle(pv_title, 'Interpreter', 'none');
print(fig, strcat(fig_path, 'Fig6CandF-PopVecCorr-history-vs-upcoming-mean_RR_LL'), '-dsvg');

%% Plot time-resolved auROC for each population
for pair_idx = 1:num_pair
    for pop_idx = 1:num_pop
        fig = figure('Position', [90 90 760 520]);
        hold on
        h_line = gobjects(numel(plot_cond_list), 1);
        for plot_idx = 1:numel(plot_cond_list)
            cond_idx = plot_cond_list(plot_idx);

            % Full population as a solid reference
            if pop_idx ~= 1
                temp_auroc = lda_results(1).cond(cond_idx).auroc(:, :, pair_idx);
                meanAUC_full = mean(temp_auroc, 2, 'omitnan');
                semAUC_full = std(temp_auroc, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(temp_auroc), 2));
                h = fill([time_bin_win'; flipud(time_bin_win')], [meanAUC_full+semAUC_full; flipud(meanAUC_full-semAUC_full)], cond_fill(cond_idx, :), 'LineStyle', 'none', 'HandleVisibility', 'off');
                set(h, 'FaceAlpha', 0.25);
            end

            temp_auroc = lda_results(pop_idx).cond(cond_idx).auroc(:, :, pair_idx);
            meanAUC = mean(temp_auroc, 2, 'omitnan');
            semAUC = std(temp_auroc, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(temp_auroc), 2));
            h = fill([time_bin_win'; flipud(time_bin_win')], [meanAUC+semAUC; flipud(meanAUC-semAUC)], cond_fill(cond_idx, :), 'LineStyle', 'none', 'HandleVisibility', 'off');
            set(h, 'FaceAlpha', 0.25);

            if pop_idx == 1
                h_line(plot_idx) = plot(time_bin_win, meanAUC, '-', 'Color', cond_line(cond_idx, :), 'LineWidth', 2.8);
            else
                h_line(plot_idx) = plot(time_bin_win, meanAUC, '--', 'Color', cond_line(cond_idx, :), 'LineWidth', 1.2);
                plot(time_bin_win, meanAUC_full, '-', 'Color', cond_line(cond_idx, :), 'LineWidth', 2.8, 'HandleVisibility', 'off');
            end
        end
        yline(0.5, 'k--', 'LineWidth', 0.8);
        xline(0, ':k', 'LineWidth', 1.3);
        xline(param.sample_time, ':k', 'LineWidth', 1.3);
        xline(param.gocue_time, ':k', 'LineWidth', 1.3);
        xlim([time_bin_win(1), time_bin_win(end)]);
        ylim([0.4, 0.8]);
        grid on
        yticks(0.4:0.1:0.8)
        xlabel('Causal window end time [s]');
        ylabel('auROC');
        title(pair_info(pair_idx).name, 'Interpreter', 'none');
        legend(h_line, cellstr(cond_name_list(plot_cond_list)), 'Location', 'southeast');
        set(gca, 'FontSize', 12);
        sgtitle(sprintf('%s | %s | %d units | R=%d', pop_info(pop_idx).label, pair_info(pair_idx).name, num_units_eff, param.num_repeat), 'Interpreter', 'none');
        print(fig, strcat(fig_path, 'Fig6CandF-', pair_info(pair_idx).short, '-timeseries-', pop_info(pop_idx).name), '-dsvg');
    end
end

%% Bar graph for the Full population (correct trials, 0.8 s vs 1.5 s within each condition)
temp_stat = stat_bar{1, 1};

fig = figure('Position', [120 120 720 540]);
hold on
Plot_bar_with_hatch(gca, temp_stat.mean(plot_cond_list, :), temp_stat.sem(plot_cond_list, :), cond_line(plot_cond_list, :));

tick = 0.018*diff(param.bar_ylim);
temp_top = temp_stat.mean(plot_cond_list, :) + temp_stat.sem(plot_cond_list, :);
y_pos = max(temp_top(:)) + 2*tick;
for plot_idx = 1:num_plot_cond
    x_pos = 2*(plot_idx-1)+[1, 2];
    Plot_sig_bracket(gca, x_pos(1), x_pos(2), y_pos, temp_stat.p_within(plot_cond_list(plot_idx)), tick, 11);
end

xticks(1:2*num_plot_cond);
xticklabels(bar_labels);
xtickangle(35);
xlim([0.3, 2*num_plot_cond+0.7]);
ylim([param.bar_ylim(1), max(param.bar_ylim(2), y_pos+3*tick)]);
yticks(param.bar_ylim(1):0.1:param.bar_ylim(2));
ylabel('auROC');
title(sprintf('%s | %s | auROC at %.2f s (plain) vs %.2f s (hatched) | %d units | R=%d', pop_info(1).label, pair_info(1).name, bar_time_actual(1), bar_time_actual(2), num_units_eff, param.num_repeat), 'Interpreter', 'none', 'FontSize', 10);
set(gca, 'FontSize', 12);
print(fig, strcat(fig_path, 'Fig6CandF-correct-bar0_8vs1_5-Full'), '-dsvg');

%% Plot time-resolved auROC for the statistics-removed population only (correct trials)
fig = figure('Position', [90 90 760 520]);
hold on
h_line = gobjects(numel(plot_cond_list), 1);
for plot_idx = 1:numel(plot_cond_list)
    cond_idx = plot_cond_list(plot_idx);
    temp_auroc = lda_results(2).cond(cond_idx).auroc(:, :, 1);
    meanAUC = mean(temp_auroc, 2, 'omitnan');
    semAUC = std(temp_auroc, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(temp_auroc), 2));
    h = fill([time_bin_win'; flipud(time_bin_win')], [meanAUC+semAUC; flipud(meanAUC-semAUC)], cond_fill(cond_idx, :), 'LineStyle', 'none', 'HandleVisibility', 'off');
    set(h, 'FaceAlpha', 0.25);
    h_line(plot_idx) = plot(time_bin_win, meanAUC, '-', 'Color', cond_line(cond_idx, :), 'LineWidth', 2.8);
end
yline(0.5, 'k--', 'LineWidth', 0.8);
xline(0, ':k', 'LineWidth', 1.3);
xline(param.sample_time, ':k', 'LineWidth', 1.3);
xline(param.gocue_time, ':k', 'LineWidth', 1.3);
xlim([time_bin_win(1), time_bin_win(end)]);
ylim([0.4, 0.8]);
grid on
yticks(0.4:0.1:0.8);
xlabel('Causal window end time [s]');
ylabel('auROC');
title(pair_info(1).name, 'Interpreter', 'none');
legend(h_line, cellstr(cond_name_list(plot_cond_list)), 'Location', 'southeast');
set(gca, 'FontSize', 12);
sgtitle(sprintf('%s | %s | %d units | R=%d', pop_info(2).label, pair_info(1).name, num_units_eff, param.num_repeat), 'Interpreter', 'none');
print(fig, strcat(fig_path, 'Fig6CandF-correct-timeseries-dropSel_Response-only-solid'), '-dsvg');

% same figure with the y-axis up to 0.7
ylim([0.4, 0.7]);
yticks(0.4:0.1:0.7);
print(fig, strcat(fig_path, 'Fig6CandF-correct-timeseries-dropSel_Response-only-solid-ylim0_4to0_7'), '-dsvg');

%% Save statistics
save(strcat(param.save_path, 'stats_', param.version, '-Fig6CandF-matched_trials_Unit', num2str(num_units_eff), '.mat'), 'stat_bar', 'rho', 'stat_pv', 'bar_time_actual', '-v7.3')

%% Functions
function Plot_bar_with_hatch(ax, mean_vals, sem_vals, colors)
% one row per condition: plain bar for the first column, hatched bar for the second column
for cond_idx = 1:size(mean_vals, 1)
    x_pos = 2*(cond_idx-1)+[1, 2];

    h1 = bar(ax, x_pos(1), mean_vals(cond_idx, 1), 'FaceColor', colors(cond_idx, :));
    errorbar(ax, x_pos(1), mean_vals(cond_idx, 1), sem_vals(cond_idx, 1), 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 12);
    h2 = bar(ax, x_pos(2), mean_vals(cond_idx, 2), 'FaceColor', colors(cond_idx, :));
    errorbar(ax, x_pos(2), mean_vals(cond_idx, 2), sem_vals(cond_idx, 2), 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 12);
    hatchfill2(h2, 'cross', 'LineWidth', 0.5);

    set(h1, 'HandleVisibility', 'off');
    set(h2, 'HandleVisibility', 'off');
end
end

function Plot_sig_bracket(ax, x1, x2, y, p, tick, font_size)
plot(ax, [x1 x1 x2 x2], [y y+tick y+tick y], '-k', 'LineWidth', 1.1, 'HandleVisibility', 'off');
text(ax, (x1+x2)/2, y+tick, Convert_pvalue_to_asterisk(p), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', font_size);
end

function asterisk = Convert_pvalue_to_asterisk(p)
if isnan(p)
    asterisk = 'n/a';
elseif p < 0.001
    asterisk = '***';
elseif p < 0.01
    asterisk = '**';
elseif p < 0.05
    asterisk = '*';
else
    asterisk = 'n.s.';
end
end
