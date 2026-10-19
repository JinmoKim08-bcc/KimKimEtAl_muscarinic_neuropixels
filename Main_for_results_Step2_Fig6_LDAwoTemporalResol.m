clc
clear
close all

%% Hyperparameters to adjust
param = struct();

% add path NPY library for reading Phy
npy_path = "utils\npy-matlab\npy-matlab\";
addpath(npy_path)

save_path = 'C:\Users\sirom\Documents\KimKimEtAl_muscarnic_neuropixels\';
parts_path = '\fig6\';

if ~isfolder(strcat(save_path, parts_path))
    mkdir(strcat(save_path, parts_path));
end

%% Hyperparameters that most do not need to be adjusted
param.fs = 30000; % [Hz]; sampling rate
param.bin_time = 0.001; % [sec]; temporal binning size to calculate firing rate
param.smoothing_pts = 100; % [idx]; smoothing data points
param.smoothing_window = ones(1, param.smoothing_pts) / param.smoothing_pts;
param.boxcar.a = 1;
param.boxcar.b = ones(1, param.smoothing_pts);
param.ITI_time = 2; % [sec]; Inter-trial interval duration
param.sample_time = 0.1; % [sec]; timing given a sample
param.delay_time = 2; % [sec]
param.gocue_time = param.sample_time + param.delay_time; % [sec]
param.response_time = 0.5; % [sec]
param.t_max = param.sample_time+ param.delay_time + param.response_time; % [sec]; maximal end time to analyze
param.t_max_truncated = param.sample_time+ param.delay_time + param.response_time; % [sec];
param.epochs.stimulus = [0, param.sample_time];
param.epochs.choice = [param.gocue_time-0.6, param.gocue_time];
param.epochs.action = [param.gocue_time, param.gocue_time+0.5];
param.epochs.outcome = [param.gocue_time+2, param.gocue_time+3];
param.eps = 1e-1;
param.time_epochs = [param.sample_time, (param.sample_time+param.delay_time), (param.sample_time+param.delay_time+param.response_time)];
% parameters for elimination
param.poor_fr = 0.5; % [Hz]


% parameters for coding direction analysis
param.fs_bin = 1/param.bin_time; % [Hz]
param.delay_dur = 0.6;
param.response_dur = 0.4;

% postulated probe channel locations for each region
param.loc_name = ["DP"; "IL"; "PrL"; "ACC"; "Motor"]; % Region names
param.prob_loc_list = [0, 14; % DP
    15, 66; % IL
    67, 162; % PrL
    163, 228; % ACC
    229, 270]; % Motor

trial_types = ["RR"; "LL"; "RL"; "LR"];
cor_trial_types = ["RR"; "LL"];

param.save_path = save_path;

param.num_iter = 100; % number of shuffling

param.version = "v16fd";

param.do_norm = 4; % 0: do not; 1: z-score; 2: normalized by baseline mean == x / mu; 3: mean subtraction; 4: z-score without param.eps
param.num_bootstrap = 100;

%% Read calculated firing vectors
cond_name_list= ["Before"; "Atropine"; "Saline"];

num_sess_list = [5, 4, 3];

max_pts = param.t_max_truncated*1/param.bin_time; %
ITI_dur = (1:param.ITI_time*param.fs)'; % 2 [sec]
ITI_dur_bin = numel(Calculate_bins_v2(ITI_dur, param.fs, param.bin_time));

load(strcat(param.save_path, 'results_', param.version, '.mat'));

%% 
clc
close all

fpath = strcat(param.save_path, '\', param.version, '\');

sep_trials = ["cor"; "CR"; "CL"; "CR-CL"; "incor"; "ER"; "EL"; "ER-EL"; "CR-ER"; "CL-EL"; "Cor+Incor"];

time_bin = -param.ITI_time:param.bin_time:param.t_max_truncated-param.bin_time;
% time_bin = (0:param.bin_time:param.t_max-param.bin_time)';
z_score_bounds = [0, 3];
time_bounds = [-1, time_bin(end)];

stats = struct();

%% LDA analysis w/o/ time resolution - Single-Single (v0d) - 100 units
time_bin = -param.ITI_time:param.bin_time:param.t_max_truncated-param.bin_time;

% parameters for pseudo population
param.num_repeat = 100;
param.num_data = 100;
param.num_chosen_units = 100;
param.ratio_train = 0.8;
param.ratio_test = 1 - param.ratio_train;

lda_results = struct();

for cond_idx = 1:numel(cond_name_list)

    datasets = struct();
    trials = struct();
    cv_opts = struct();

    % Sample session
    rng(cond_idx, 'twister');
    num_total_sess = num_sess_list(cond_idx);
    sess_sampled = {};
    for repeat_idx = 1:param.num_repeat
        sess_list = randsample(num_total_sess, num_total_sess, true);      % replace true
        sess_sampled{repeat_idx} = sess_list;
    end

    % Sample trials
    for repeat_idx = 1:param.num_repeat
        sess_list = sess_sampled{repeat_idx};
        for temp_sess_idx = 1:numel(sess_list)
            sess_idx = sess_list(temp_sess_idx);

            for trial_idx = 1:numel(trial_types)
                selected_trial_type = trial_types(trial_idx);

                % Read sizes only
                vals = results(cond_idx).sess(sess_idx).(selected_trial_type).allTrials.vals;
                num_trials = size(vals, 3);
                num_units  = size(vals, 1);

                % Number of pseudo-population samples
                nTrain = floor(param.num_data * param.ratio_train);
                nTest  = max(1, floor(param.num_data * (1 - param.ratio_train)));  % avoid 0

                fprintf('Sess: [%d/%d], Trial: [%d/%d], Repeat: [%d/%d]\n', ...
                    sess_idx, num_sess_list(cond_idx), trial_idx, numel(trial_types), repeat_idx, param.num_repeat);

                rng(repeat_idx); % per-repeat seed (adjust if separate seeds per session/type are needed)

                % --- (Key) split the trial pool so that train/test never overlap ---
                if num_trials >= 2
                    perm = randperm(num_trials);
                    % Pool size: keep at least 1 trial in each
                    nPoolTrain = max(1, min(num_trials-1, round(num_trials * param.ratio_train)));
                    train_pool = perm(1:nPoolTrain);
                    test_pool  = perm(nPoolTrain+1:end);
                else
                    % Cannot split with only 1 trial -> both set to 1 (unavoidable)
                    train_pool = 1;
                    test_pool  = 1;
                end

                % --- Generate indices by sampling only within each pool (with replacement) ---
                % Size: [num samples x num units]
                if ~isempty(train_pool)
                    train_indices_all = uint32( ...
                        train_pool( randi(numel(train_pool), nTrain, num_units) ) );
                else
                    train_indices_all = ones(nTrain, num_units, 'uint32'); % safeguard for the abnormal case
                end

                if ~isempty(test_pool)
                    test_indices_all = uint32( ...
                        test_pool( randi(numel(test_pool), nTest, num_units) ) );
                else
                    test_indices_all = ones(nTest, num_units, 'uint32'); % safeguard for the abnormal case
                end

                % Store (rows: data, columns: units)
                trials(temp_sess_idx).(selected_trial_type).repeat(repeat_idx).trainSets = train_indices_all;
                trials(temp_sess_idx).(selected_trial_type).repeat(repeat_idx).testSets  = test_indices_all;
            end
        end
    end

    % Sample units
    for repeat_idx = 1:param.num_repeat
        sess_list = sess_sampled{repeat_idx};

        per_sess = floor(param.num_chosen_units / numel(sess_list));

        num_units_per_sess = ones(numel(sess_list), 1).*per_sess;
        remained_units = param.num_chosen_units - sum(num_units_per_sess);
        num_units_per_sess(1) = num_units_per_sess(1) + remained_units;

        for temp_sess_idx = 1:numel(sess_list)
            sess_idx = sess_list(temp_sess_idx);
            vals = results(cond_idx).sess(sess_idx).RR.allTrials.vals;
            num_units  = size(vals, 1);

            unit_list = randsample(num_units, num_units_per_sess(temp_sess_idx), false);      % replace true

            trials(temp_sess_idx).units.repeat(repeat_idx).units = unit_list;
        end
    end

    num_units_cumsum2 = cumsum(num_units_per_sess);
    num_units_cumsum1 = [1; num_units_cumsum2(1:end-1)+1];
    num_units_cumsum = [num_units_cumsum1, num_units_cumsum2];

    % Make a pseudo population
    % --- Precompute fixed indices/sizes ---
    nTrain = floor(param.num_data * param.ratio_train);
    nTest  = floor(param.num_data * param.ratio_test);
    twin   = 2001:numel(time_bin);          % time window to average over
    Utot   = param.num_chosen_units;       % total number of units

    for repeat_idx = 1:param.num_repeat
        sess_list = sess_sampled{repeat_idx};

        for trial_idx = 1:numel(trial_types)
            selected_trial_type = trial_types(trial_idx);

            % Output buffers: final features are already time-averaged (1 value/unit), so only 2D is needed
            train_data = zeros(nTrain, Utot, 'single');
            test_data  = zeros(nTest , Utot, 'single');

            % ---------- Do "preprocessing + time averaging" only once per session ----------
            % For each session s: vals [U x T x R] -> smooth -> baseline subtract -> time-avg
            % Result: avg_by_trial{s}: [U_s x R_s] (units x trials)
            avg_by_trial = cell(numel(sess_list), 1);

            for temp_sess_idx = 1:numel(sess_list)
                sess_idx = sess_list(temp_sess_idx);
                vals = results(cond_idx).sess(sess_idx).(selected_trial_type).allTrials.vals;   % [U x T x R]
                [U, T, R] = size(vals);

                selected_unit_indices = trials(temp_sess_idx).units.repeat(repeat_idx).units;
                vals = vals(selected_unit_indices, :, :);

                % Smoothing along the time axis (dim=2)
                vals = movmean(vals, param.smoothing_pts, 2);

                % Subtract the baseline (1:ITI_dur_bin) mean: implicit expansion, no bsxfun
                mu0 = mean(vals(:, 1:ITI_dur_bin, :), 2);     % [U x 1 x R]
                vals = vals - mu0;                             % [U x T x R]

                % Only the mean over 2001:end is used in the end
                vavg = mean(vals(:, twin, :), 2);              % [U x 1 x R]
                avg_by_trial{temp_sess_idx} = squeeze(vavg);  % [U x R]
            end

            % ---------- Build the training set (only selects trials per unit) ----------
            % (Inside the repeat loop) build the training set after per-session preprocessing
            for data_idx = 1:nTrain
                if mod(data_idx, 10) == 0
                    fprintf('Repeat: [%d/%d], Trial: [%d/%d], Data(train): [%d/%d]\n', ...
                        repeat_idx, param.num_repeat, trial_idx, numel(trial_types), data_idx, nTrain);
                end

                row = zeros(1, Utot, 'single');
                for temp_sess_idx = 1:numel(sess_list)
                    lb = num_units_cumsum(temp_sess_idx, 1);
                    ub = num_units_cumsum(temp_sess_idx, 2);

                    A = avg_by_trial{temp_sess_idx};      % [U_s x R_s]
                    % Selected unit indices (session-local unit numbers)
                    selected_unit_indices = trials(temp_sess_idx).units.repeat(repeat_idx).units;
                    U = numel(selected_unit_indices);     % number of selected units = number of rows of A

                    % ★ Subset trainSets so that only the 'selected unit' columns are taken ★
                    idx_per_unit = trials(temp_sess_idx).(selected_trial_type) ...
                        .repeat(repeat_idx).trainSets(data_idx, selected_unit_indices).';   % [U x 1]

                    % Now (1:U).' and idx_per_unit have the same length
                    lin = sub2ind(size(A), (1:U).', idx_per_unit);
                    row(lb:ub) = A(lin);
                end
                train_data(data_idx, :) = row;
            end

            % Apply the same change to the test set
            for data_idx = 1:nTest
                row = zeros(1, Utot, 'single');
                for temp_sess_idx = 1:numel(sess_list)
                    lb = num_units_cumsum(temp_sess_idx, 1);
                    ub = num_units_cumsum(temp_sess_idx, 2);

                    A = avg_by_trial{temp_sess_idx};
                    selected_unit_indices = trials(temp_sess_idx).units.repeat(repeat_idx).units;
                    U = numel(selected_unit_indices);

                    idx_per_unit = trials(temp_sess_idx).(selected_trial_type) ...
                        .repeat(repeat_idx).testSets(data_idx, selected_unit_indices).';    % [U x 1]

                    lin = sub2ind(size(A), (1:U).', idx_per_unit);
                    row(lb:ub) = A(lin);
                end
                test_data(data_idx, :) = row;
            end

            % Store
            datasets(repeat_idx).(selected_trial_type).train_data = train_data;
            datasets(repeat_idx).(selected_trial_type).test_data  = test_data;
        end
    end


    fpr_list = {};
    tpr_list = {};
    for repeat_idx = 1:param.num_repeat
        %%% Correct trials
        selected_trial_type = "RR";
        X_train_R = datasets(repeat_idx).(selected_trial_type).train_data; % N x T
        X_test_R = datasets(repeat_idx).(selected_trial_type).test_data; % N x T
        selected_trial_type = "LL";
        X_train_L = datasets(repeat_idx).(selected_trial_type).train_data; % N x T
        X_test_L = datasets(repeat_idx).(selected_trial_type).test_data; % N x T
        [R, N] = size(X_train_R);

        X_train = [X_train_R; X_train_L]; % 2R x N
        Y_train = [repmat("Right", R, 1) ; ...
                   repmat("Left" , R, 1)]; % 2R x 1

        ldaMdl = fitcdiscr(X_train, Y_train, 'DiscrimType','pseudoLinear', 'Prior','empirical');
        % ldaMdl = fitcdiscr(X_train, Y_train, 'DiscrimType','linear', 'Gamma', 0.1, 'Delta', 0, 'Prior', 'empirical');

        % Extract the discriminant vector between class 1 vs 2
        w = ldaMdl.Coeffs(1,2).Linear;   % [N×1] discriminant vector
        b = ldaMdl.Coeffs(1,2).Const;   % scalar intercept

        [R, N] = size(X_test_R);

        X_test = [X_test_R; X_test_L]; % 2T x N
        Y_test = [repmat("Right", R, 1) ; ...
                  repmat("Left" , R, 1)]; % 2T x 1

        proj_test = X_test*w + b;

        % proj_cor_R = [proj_cor_R, proj_test(Y_test=="Right")];
        % proj_cor_L = [proj_cor_L, proj_test(Y_test=="Left" )];

        %%% metrics
        [Y_hat, score] = predict(ldaMdl, X_test);
        accuracy  = mean(Y_hat == Y_test);
        posLabel    = {'Right'};   % e.g., treat the first class as positive
        [isPos, posIdx] = ismember(posLabel, ldaMdl.ClassNames);
        [fpr,tpr,~,auROC] = perfcurve(Y_test, score(:,posIdx), posLabel);
        % acc_cor = [acc_cor; accuracy];
        % auroc_cor = [auroc_cor; auROC];

        fpr_list_cor{repeat_idx} = fpr;
        tpr_list_cor{repeat_idx} = tpr;

        %%% Incorrect trials
        selected_trial_type = "LR";
        X_train_R = datasets(repeat_idx).(selected_trial_type).train_data; % N x T
        X_test_R = datasets(repeat_idx).(selected_trial_type).test_data; % N x T
        selected_trial_type = "RL";
        X_train_L = datasets(repeat_idx).(selected_trial_type).train_data; % N x T
        X_test_L = datasets(repeat_idx).(selected_trial_type).test_data; % N x T
        [R, N] = size(X_train_R);

        X_train = [X_train_R; X_train_L]; % 2R x N
        Y_train = [repmat("Right", R, 1) ; ...
                   repmat("Left" , R, 1)]; % 2R x 1

        ldaMdl = fitcdiscr(X_train, Y_train, 'DiscrimType','pseudoLinear', 'Prior','empirical');
        % ldaMdl = fitcdiscr(X_train, Y_train, 'DiscrimType','linear', 'Gamma', 0.1, 'Delta', 0, 'Prior', 'empirical');

        % Extract the discriminant vector between class 1 vs 2
        w = ldaMdl.Coeffs(1,2).Linear;   % [N×1] discriminant vector
        b = ldaMdl.Coeffs(1,2).Const;   % scalar intercept

        [R, N] = size(X_test_R);

        X_test = [X_test_R; X_test_L]; % 2T x N
        Y_test = [repmat("Right", R, 1) ; ...
                  repmat("Left" , R, 1)]; % 2T x 1

        % proj_test = X_test*w + b;

        % proj_incor_R = [proj_incor_R, proj_test(Y_test=="Right")];
        % proj_incor_L = [proj_incor_L, proj_test(Y_test=="Left" )];

        %%% metrics
        [Y_hat, score] = predict(ldaMdl, X_test);
        accuracy  = mean(Y_hat == Y_test);
        posLabel    = {'Right'};   % e.g., treat the first class as positive
        [isPos, posIdx] = ismember(posLabel, ldaMdl.ClassNames);
        [fpr,tpr,~,auROC] = perfcurve(Y_test, score(:,posIdx), posLabel);
        % acc_incor = [acc_incor; accuracy];
        % auroc_incor = [auroc_incor; auROC];

        fpr_list_incor{repeat_idx} = fpr;
        tpr_list_incor{repeat_idx} = tpr;
    end

    lda_results(cond_idx).fpr_list_cor = fpr_list_cor;
    lda_results(cond_idx).tpr_list_cor = tpr_list_cor;
    lda_results(cond_idx).fpr_list_incor = fpr_list_incor;
    lda_results(cond_idx).tpr_list_incor = tpr_list_incor;
end

save(strcat(param.save_path, '\', 'selective_index_', param.version, '-0toEnd_PopUnits-single-single_Unit', num2str(param.num_chosen_units), '.mat'), 'lda_results', '-v7.3')
load(strcat(param.save_path, '\', 'selective_index_', param.version, '-0toEnd_PopUnits-single-single_Unit', num2str(param.num_chosen_units), '.mat'))

%% Plot for ROC curve (w/o/ temporal resolution)
for cond_idx = 1:numel(cond_name_list)

    fpr_list_cor = lda_results(cond_idx).fpr_list_cor;
    tpr_list_cor = lda_results(cond_idx).tpr_list_cor;
    fpr_list_incor = lda_results(cond_idx).fpr_list_incor;
    tpr_list_incor = lda_results(cond_idx).tpr_list_incor;

    %%% ROC curve for correct
    [fpr_cor, tpr_cor] = Interpolate_ROC_curve(fpr_list_cor, tpr_list_cor, param);
    meanTPR = mean(tpr_cor,   1)';
    semTPR  = std(tpr_cor, 0, 1)' ./ sqrt(param.num_repeat);

    fig = figure;
    hold on;
    plot(fpr_cor, meanTPR, 'LineWidth',2, 'color', 'g');
    h = fill([fpr_cor; flipud(fpr_cor)], [meanTPR+semTPR; flipud(meanTPR-semTPR)], 'g', 'EdgeColor','none');
    set(h, 'FaceAlpha', 0.2);

    %%% ROC curve for incorrect
    [fpr_incor, tpr_incor] = Interpolate_ROC_curve(fpr_list_incor, tpr_list_incor, param);
    meanTPR = mean(tpr_incor,   1)';
    semTPR  = std(tpr_incor, 0, 1)' ./ sqrt(param.num_repeat);

    plot(fpr_incor, meanTPR, 'LineWidth',2, 'color', '#ff8c00');
    h = fill([fpr_incor; flipud(fpr_incor)], [meanTPR+semTPR; flipud(meanTPR-semTPR)], 'r', 'EdgeColor','none');
    set(h, 'FaceColor', "#ff8c00")
    set(h, 'FaceAlpha', 0.2);

    % Diagonal (chance)
    plot([0 1],[0 1],'k--');
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    title('Average ROC with \pm SEM');
    hold off;
    xlim([0, 1])
    ylim([0, 1])
    set(gca, 'FontSize', 15)
    print(fig, strcat(param.save_path, parts_path, cond_name_list(cond_idx), '-LDA_wo_Time_resol_AllUnits_single-single-Unit', num2str(param.num_chosen_units)), '-dsvg');
end

%% Plot for ROC curve (w/o/ temporal resolution) & hold on
clc
close all

fig = figure;
hold on

selected_colors = {};
selected_colors{1,1} = [0.5, 1, 0.5];
selected_colors{1,2} = [0.5, 1, 0.5];

selected_colors{2,1} = [1.0, 0.549, 0.0];
selected_colors{2,2} = [1.0, 0.549, 0.0];

selected_colors{3,1} = [0.494, 0.184, 0.557];
selected_colors{3,2} = [0.494, 0.184, 0.557];

for cond_idx = 1:numel(cond_name_list)
    fpr_list_cor = lda_results(cond_idx).fpr_list_cor;
    tpr_list_cor = lda_results(cond_idx).tpr_list_cor;

    [fpr_cor, tpr_cor] = Interpolate_ROC_curve(fpr_list_cor, tpr_list_cor, param);

    % Mean/SD ignoring NaNs
    meanTPR = mean(tpr_cor, 1, 'omitnan')';        % N×1
    semTPR  = std(tpr_cor, 0, 1, 'omitnan')' ./ sqrt(param.num_repeat);

    % Convert to a column vector
    fpr_cor = fpr_cor(:);

    % Line
    plot(fpr_cor, meanTPR, 'LineWidth', 2, 'Color', selected_colors{cond_idx,1});

    % Band (filled area): 2N×1 closed polygon
    x_fill = [fpr_cor; flipud(fpr_cor)];
    y_fill = [meanTPR + semTPR; flipud(meanTPR - semTPR)];
    h = fill(x_fill, y_fill, selected_colors{cond_idx,2}, 'EdgeColor','none');
    set(h, 'FaceAlpha', 0.2);
end

% Reference line, labels, and axis limits are set only once, outside the loop
plot([0 1],[0 1],'k--');              % Chance
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('Average ROC with \pm SEM');
xlim([0, 1]); ylim([0, 1]);
set(gca, 'FontSize', 15);
hold off

print(fig, strcat(param.save_path, parts_path, 'HoldON', '-LDA_wo_Time_resol_AllUnits_single-single-Unit', num2str(param.num_chosen_units)), '-dsvg');

%% auROC bar plot (w/o/ temporal resolution): bootstrap mean +/- SEM
% Correct trials, matching the overlaid ROC curves above.
% Each original ROC contributes one auROC (100 per condition by default).
% Integrate BEFORE interpolation: removing duplicate FPRs can change the area.
auroc_cond_names = ["Before"; "Saline"; "Atropine"];
[has_condition, auroc_cond_indices] = ismember(lower(auroc_cond_names), ...
    lower(string(cond_name_list(:))));
assert(all(has_condition), 'Before, Saline, and Atropine results are required.');

num_auroc_conditions = numel(auroc_cond_names);
auroc_values = nan(param.num_repeat, num_auroc_conditions);
for plot_idx = 1:num_auroc_conditions
    cond_idx = auroc_cond_indices(plot_idx);
    fpr_bootstrap = lda_results(cond_idx).fpr_list_cor;
    tpr_bootstrap = lda_results(cond_idx).tpr_list_cor;
    assert(numel(fpr_bootstrap) == param.num_repeat && ...
        numel(tpr_bootstrap) == param.num_repeat, ...
        'Expected %d bootstrap ROC curves for %s.', ...
        param.num_repeat, char(auroc_cond_names(plot_idx)));

    for repeat_idx = 1:param.num_repeat
        fpr_boot = fpr_bootstrap{repeat_idx}(:);
        tpr_boot = tpr_bootstrap{repeat_idx}(:);
        assert(numel(fpr_boot) >= 2 && numel(fpr_boot) == numel(tpr_boot) && ...
            all(isfinite(fpr_boot)) && all(isfinite(tpr_boot)) && ...
            all(fpr_boot >= 0 & fpr_boot <= 1) && ...
            all(tpr_boot >= 0 & tpr_boot <= 1) && ...
            all(diff(fpr_boot) >= 0) && all(diff(tpr_boot) >= 0) && ...
            fpr_boot(1) == 0 && fpr_boot(end) == 1, ...
            'Invalid ROC for %s, bootstrap %d.', ...
            char(auroc_cond_names(plot_idx)), repeat_idx);
        auroc_values(repeat_idx, plot_idx) = trapz(fpr_boot, tpr_boot);
    end
    lda_results(cond_idx).auroc_cor = auroc_values(:, plot_idx);
end

auroc_n = sum(isfinite(auroc_values), 1);
assert(all(auroc_n >= 2), 'At least two bootstrap auROCs per condition are required.');
auroc_mean = mean(auroc_values, 1);
auroc_sem = std(auroc_values, 0, 1) ./ sqrt(auroc_n);
% This SEM describes the mean across bootstrap replicates, not biological SEM.
% Bootstrap replicates are resamples, not independent animals/sessions.
auroc_summary = table(auroc_cond_names, auroc_n(:), auroc_mean(:), auroc_sem(:), ...
    'VariableNames', {'Condition', 'N_bootstrap', 'Mean_auROC', 'SEM_auROC'});

% Two-sided Mann-Whitney U (Wilcoxon rank-sum), corrected across all 3 pairs.
auroc_pairs = [1 2; 1 3; 2 3];
num_auroc_comparisons = size(auroc_pairs, 1);
auroc_p_raw = nan(num_auroc_comparisons, 1);
auroc_u = nan(num_auroc_comparisons, 1);
for pair_idx = 1:num_auroc_comparisons
    idx1 = auroc_pairs(pair_idx, 1);
    idx2 = auroc_pairs(pair_idx, 2);
    pair_values = auroc_values(:, [idx1 idx2]);
    if all(pair_values(:) == pair_values(1))
        % ranksum's normal approximation has zero variance when all ranks tie.
        auroc_p_raw(pair_idx) = 1;
        auroc_u(pair_idx) = auroc_n(idx1)*auroc_n(idx2)/2;
        continue
    end
    [auroc_p_raw(pair_idx), ~, rank_stats] = ranksum( ...
        auroc_values(:, idx1), auroc_values(:, idx2), 'tail', 'both');
    auroc_u(pair_idx) = rank_stats.ranksum - auroc_n(idx1)*(auroc_n(idx1)+1)/2;
end
assert(all(isfinite(auroc_p_raw)), 'Mann-Whitney U test returned a nonfinite p value.');
auroc_p_bonferroni = min(auroc_p_raw * num_auroc_comparisons, 1);
auroc_significance = repmat("ns", num_auroc_comparisons, 1);
auroc_significance(auroc_p_bonferroni < 0.05) = "*";
auroc_significance(auroc_p_bonferroni < 0.01) = "**";
auroc_significance(auroc_p_bonferroni < 0.001) = "***";
auroc_stats = table(auroc_cond_names(auroc_pairs(:, 1)), ...
    auroc_cond_names(auroc_pairs(:, 2)), auroc_u, auroc_p_raw, ...
    auroc_p_bonferroni, auroc_significance, ...
    'VariableNames', {'Condition1', 'Condition2', 'U_condition1', ...
    'P_raw', 'P_Bonferroni', 'Significance'});
disp(auroc_summary);
disp(auroc_stats);

fig_auroc = figure('Color', 'w');
hold on
% Preserve the condition colors from the ROC overlay: green, purple, orange.
auroc_colors = [0.5, 1, 0.5; 0.494, 0.184, 0.557; 1.0, 0.549, 0.0];
% Use one opaque, single-face patch per bar for individual vector editing.
% Keep the original bar width and zero baseline (clipped by the y-axis limits).
auroc_bar_width = 0.65;
bar_auroc = gobjects(num_auroc_conditions, 1);
for plot_idx = 1:num_auroc_conditions
    bar_x = plot_idx + [-1, 1, 1, -1]*auroc_bar_width/2;
    bar_y = [0, 0, auroc_mean(plot_idx), auroc_mean(plot_idx)];
    bar_auroc(plot_idx) = patch(bar_x, bar_y, auroc_colors(plot_idx, :), ...
        'FaceAlpha', 1, 'EdgeColor', 'k', 'LineWidth', 0.5);
end
errorbar(1:num_auroc_conditions, auroc_mean, auroc_sem, ...
    'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 12);

% Draw short brackets first and the Before-Atropine bracket above them.
auroc_bracket_order = [1 3 2];
% Keep all three brackets and labels inside the fixed [0.4, 1] y-axis.
auroc_bracket_step = 0.05;
auroc_bracket_top = 0.94;
auroc_bracket_base = min(max(auroc_mean + auroc_sem) + 0.03, ...
    auroc_bracket_top - (num_auroc_comparisons-1)*auroc_bracket_step);
for bracket_idx = 1:num_auroc_comparisons
    pair_idx = auroc_bracket_order(bracket_idx);
    x_pair = auroc_pairs(pair_idx, :);
    bracket_y = auroc_bracket_base + (bracket_idx-1)*auroc_bracket_step;
    plot([x_pair(1), x_pair(1), x_pair(2), x_pair(2)], ...
        [bracket_y-0.015, bracket_y, bracket_y, bracket_y-0.015], ...
        'k-', 'LineWidth', 1.2);
    text(mean(x_pair), bracket_y+0.012, auroc_significance(pair_idx), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 13);
end
set(gca, 'XTick', 1:num_auroc_conditions, ...
    'XTickLabel', cellstr(auroc_cond_names), 'FontSize', 15);
xlim([0.5, num_auroc_conditions+0.5]);
ylim([0.4, 1]);
yticks(0.4:0.1:1);
ylabel('auROC');
title({'Correct-trial auROC (mean \pm SEM)'});
box off
hold off

auroc_output_base = strcat(param.save_path, parts_path, ...
    'auROC-Bar-LDA_wo_Time_resol_AllUnits_single-single-Unit', num2str(param.num_chosen_units));
% Force vector output so Illustrator receives editable paths, not an image.
print(fig_auroc, auroc_output_base, '-dsvg', '-vector');
save(strcat(auroc_output_base, '.mat'), 'auroc_values', 'auroc_cond_names', ...
    'auroc_summary', 'auroc_stats');
writetable(auroc_stats, strcat(auroc_output_base, '-stats.csv'));

