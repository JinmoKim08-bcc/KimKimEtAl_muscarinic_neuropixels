clc
clear
close all

%% Hyperparameters to adjust
npy_path = "utils\npy-matlab\npy-matlab\";
addpath(npy_path)

param = struct();

save_path = 'C:\Users\sirom\Documents\KimKimEtAl_muscarnic_neuropixels\';
parts_path = '\fig6\';


%% Hyperparameters that most do not need to be adjusted
param.bin_time = 0.001; % [sec]; temporal binning size to calculate firing rate
param.ITI_time = 2.0; % [sec]; Inter-trial interval duration
param.sample_time = 0.1; % [sec]; timing given a sample
param.delay_time = 2.0; % [sec]
param.gocue_time = param.sample_time+param.delay_time; % [sec]
param.response_time = 0.5; % [sec]
param.t_max = param.gocue_time+param.response_time; % [sec]; maximal end time to analyze

param.smoothing_pts = 100; % [idx]; causal smoothing data points
param.zscore_eps = 1e-6; % added to the baseline SD for z-score

% parameters for histograms
param.hist_bin_widths = 0.02*ones(5, 1); % bin width for each preprocessing
param.abs_diff_xlim = [0, 2]; % x range of the cropped |RR-LL| figures
param.hist_display_style = 'bar';
param.hist_face_alpha = 0.24;
param.hist_edge_alpha = 0.65;

% epochs
param.epoch_names = ["Stimulus"; "DelayEarly"; "DelayLate"; "Response"];
param.epoch_windows = [0.0, param.sample_time; % [sec]
    param.sample_time, param.sample_time+1.0;
    param.sample_time+1, param.gocue_time;
    param.gocue_time, param.t_max];

% preprocessing methods
param.prep_names = ["CommonBaseline"; % epoch mean of each trial - common RR+LL baseline
    "B2PCausalBaselineSub"; % trial average + causal smoothing + baseline subtraction
    "B2PCausalZScore"; % trial average + causal smoothing + baseline z-score
    "TrialAvgZScore"; % trial average + baseline z-score
    "TrialAvgBaselineSub"]; % trial average + baseline subtraction
param.prep_labels = ["Common RR+LL baseline subtraction";
    "Fig5B2P trial-average + causal smoothing + baseline subtraction";
    "Fig5B2P trial-average + causal smoothing + baseline z-score";
    "Trial-average + baseline z-score (no smoothing)";
    "Trial-average + baseline subtraction (no smoothing)"];
param.prep_units = ["Hz"; "Hz"; "z"; "z"; "Hz"];
param.prep_file_tags = ["CommonBaseline"; "B2PBaselineSub"; "B2PZScore"; "TrialAvgZScore"; "TrialAvgBaselineSub"];

cond_name_list = ["Before"; "Atropine"; "Saline"];

param.cond_name_list = cond_name_list;
param.cond_colors = [0.00, 0.65, 0.10; % Before
    1.00, 0.55, 0.00; % Atropine
    0.49, 0.18, 0.56]; % Saline
param.cond_pairs = [1, 2; 1, 3; 2, 3];
param.bar_cond_order = [1, 3, 2]; % bar order: Before / Saline / Atropine

param.save_path = save_path;
param.version = "v16fd";

%% Read calculated firing vectors
load(strcat(param.save_path, 'results_', param.version, '.mat'), 'results');

num_epochs = numel(param.epoch_names);
num_prep = numel(param.prep_names);

time_bin = -param.ITI_time:param.bin_time:(param.t_max-param.bin_time);
num_pts = numel(time_bin);
baseline_bins = find(time_bin >= -param.ITI_time & time_bin < 0);
epoch_bins = cell(num_epochs, 1);
for epoch_idx = 1:num_epochs
    epoch_bins{epoch_idx} = find(time_bin >= param.epoch_windows(epoch_idx, 1) & time_bin < param.epoch_windows(epoch_idx, 2));
end

%% RR, LL and |RR-LL| of each unit for every preprocessing and epoch
baseline_table = table();
value_table = table();

for cond_idx = 1:numel(cond_name_list)
    num_sess = numel(results(cond_idx).sess);
    cond_unit_offset = 0;

    for sess_idx = 1:num_sess
        vals_RR = results(cond_idx).sess(sess_idx).RR.allTrials.vals;
        vals_LL = results(cond_idx).sess(sess_idx).LL.allTrials.vals;
        vals_RR = vals_RR(:, 1:num_pts, :);
        vals_LL = vals_LL(:, 1:num_pts, :);
        num_units = size(vals_RR, 1);
        num_trials_RR = size(vals_RR, 3);
        num_trials_LL = size(vals_LL, 3);

        fprintf('Cond: [%d/%d], Sess: [%d/%d], Units: %d, Trials (RR/LL): %d/%d\n', cond_idx, numel(cond_name_list), sess_idx, num_sess, num_units, num_trials_RR, num_trials_LL);

        %%% Common baseline over all RR and LL trials (for CommonBaseline)
        temp_baseline = double(vals_RR(:, baseline_bins, :));
        is_finite_RR = isfinite(temp_baseline);
        temp_baseline(~is_finite_RR) = 0;
        baseline_sum_RR = sum(temp_baseline, [2, 3]);

        temp_baseline = double(vals_LL(:, baseline_bins, :));
        is_finite_LL = isfinite(temp_baseline);
        temp_baseline(~is_finite_LL) = 0;
        baseline_sum_LL = sum(temp_baseline, [2, 3]);

        common_baseline_mean = (baseline_sum_RR+baseline_sum_LL)./(sum(is_finite_RR, [2, 3])+sum(is_finite_LL, [2, 3]));

        %%% Trial-averaged PSTH (summed in single precision)
        psth_RR = zeros(num_units, num_pts, 'single');
        for trial_idx = 1:num_trials_RR
            psth_RR = psth_RR + single(vals_RR(:, :, trial_idx));
        end
        psth_RR = double(psth_RR)/num_trials_RR;

        psth_LL = zeros(num_units, num_pts, 'single');
        for trial_idx = 1:num_trials_LL
            psth_LL = psth_LL + single(vals_LL(:, :, trial_idx));
        end
        psth_LL = double(psth_LL)/num_trials_LL;

        % causal smoothing (trailing window)
        smoothed_RR = movmean(psth_RR, [param.smoothing_pts-1, 0], 2, 'Endpoints', 'shrink');
        smoothed_LL = movmean(psth_LL, [param.smoothing_pts-1, 0], 2, 'Endpoints', 'shrink');

        %%% Baseline mean, SD and z-score (RR and LL separately)
        mu_smoothed_RR = mean(smoothed_RR(:, baseline_bins), 2, 'omitnan');
        sd_smoothed_RR = std(smoothed_RR(:, baseline_bins), 0, 2, 'omitnan');
        denom_smoothed_RR = sd_smoothed_RR + param.zscore_eps;
        z_smoothed_RR = (smoothed_RR - mu_smoothed_RR)./denom_smoothed_RR;

        mu_smoothed_LL = mean(smoothed_LL(:, baseline_bins), 2, 'omitnan');
        sd_smoothed_LL = std(smoothed_LL(:, baseline_bins), 0, 2, 'omitnan');
        denom_smoothed_LL = sd_smoothed_LL + param.zscore_eps;
        z_smoothed_LL = (smoothed_LL - mu_smoothed_LL)./denom_smoothed_LL;

        mu_unsmoothed_RR = mean(psth_RR(:, baseline_bins), 2, 'omitnan');
        sd_unsmoothed_RR = std(psth_RR(:, baseline_bins), 0, 2, 'omitnan');
        denom_unsmoothed_RR = sd_unsmoothed_RR + param.zscore_eps;
        z_unsmoothed_RR = (psth_RR - mu_unsmoothed_RR)./denom_unsmoothed_RR;

        mu_unsmoothed_LL = mean(psth_LL(:, baseline_bins), 2, 'omitnan');
        sd_unsmoothed_LL = std(psth_LL(:, baseline_bins), 0, 2, 'omitnan');
        denom_unsmoothed_LL = sd_unsmoothed_LL + param.zscore_eps;
        z_unsmoothed_LL = (psth_LL - mu_unsmoothed_LL)./denom_unsmoothed_LL;

        % baseline subtraction
        sub_smoothed_RR = smoothed_RR - mu_smoothed_RR;
        sub_smoothed_LL = smoothed_LL - mu_smoothed_LL;
        sub_unsmoothed_RR = psth_RR - mu_unsmoothed_RR;
        sub_unsmoothed_LL = psth_LL - mu_unsmoothed_LL;

        %%% Baseline parameters
        sess_id = string(sprintf('%s_S%02d', cond_name_list(cond_idx), sess_idx));
        cond_unit_list = cond_unit_offset + (1:num_units)';

        temp_table = table(repmat(cond_name_list(cond_idx), num_units, 1), repmat(cond_idx, num_units, 1), repmat(sess_id, num_units, 1), repmat(sess_idx, num_units, 1), ...
            (1:num_units)', cond_unit_list, repmat(num_trials_RR, num_units, 1), repmat(num_trials_LL, num_units, 1), common_baseline_mean, ...
            mu_smoothed_RR, mu_smoothed_LL, sd_smoothed_RR, sd_smoothed_LL, denom_smoothed_RR, denom_smoothed_LL, ...
            mu_unsmoothed_RR, mu_unsmoothed_LL, sd_unsmoothed_RR, sd_unsmoothed_LL, denom_unsmoothed_RR, denom_unsmoothed_LL, ...
            'VariableNames', {'Condition', 'ConditionIndex', 'SessionID', 'Session', 'SessionUnit', 'ConditionUnit', 'NTrialsRR', 'NTrialsLL', 'CommonBaselineMeanFRHz', ...
            'SmoothedBaselineMeanRRHz', 'SmoothedBaselineMeanLLHz', 'SmoothedBaselineSDRRHz', 'SmoothedBaselineSDLLHz', 'SmoothedZDenominatorRRHz', 'SmoothedZDenominatorLLHz', ...
            'UnsmoothedBaselineMeanRRHz', 'UnsmoothedBaselineMeanLLHz', 'UnsmoothedBaselineSDRRHz', 'UnsmoothedBaselineSDLLHz', 'UnsmoothedZDenominatorRRHz', 'UnsmoothedZDenominatorLLHz'});
        baseline_table = [baseline_table; temp_table];

        %%% Epoch mean for each preprocessing
        for epoch_idx = 1:num_epochs
            selected_bins = epoch_bins{epoch_idx};

            % epoch mean within each trial (for CommonBaseline)
            epoch_mean_RR = reshape(mean(double(vals_RR(:, selected_bins, :)), 2, 'omitnan'), num_units, num_trials_RR);
            epoch_mean_LL = reshape(mean(double(vals_LL(:, selected_bins, :)), 2, 'omitnan'), num_units, num_trials_LL);

            prep_vals_RR = zeros(num_units, num_prep);
            prep_vals_RR(:, 1) = mean(epoch_mean_RR-common_baseline_mean, 2, 'omitnan'); % CommonBaseline
            prep_vals_RR(:, 2) = mean(sub_smoothed_RR(:, selected_bins), 2, 'omitnan'); % B2PCausalBaselineSub
            prep_vals_RR(:, 3) = mean(z_smoothed_RR(:, selected_bins), 2, 'omitnan'); % B2PCausalZScore
            prep_vals_RR(:, 4) = mean(z_unsmoothed_RR(:, selected_bins), 2, 'omitnan'); % TrialAvgZScore
            prep_vals_RR(:, 5) = mean(sub_unsmoothed_RR(:, selected_bins), 2, 'omitnan'); % TrialAvgBaselineSub

            prep_vals_LL = zeros(num_units, num_prep);
            prep_vals_LL(:, 1) = mean(epoch_mean_LL-common_baseline_mean, 2, 'omitnan');
            prep_vals_LL(:, 2) = mean(sub_smoothed_LL(:, selected_bins), 2, 'omitnan');
            prep_vals_LL(:, 3) = mean(z_smoothed_LL(:, selected_bins), 2, 'omitnan');
            prep_vals_LL(:, 4) = mean(z_unsmoothed_LL(:, selected_bins), 2, 'omitnan');
            prep_vals_LL(:, 5) = mean(sub_unsmoothed_LL(:, selected_bins), 2, 'omitnan');

            for prep_idx = 1:num_prep
                temp_table = table(repmat(param.prep_names(prep_idx), num_units, 1), repmat(param.prep_labels(prep_idx), num_units, 1), repmat(param.prep_units(prep_idx), num_units, 1), ...
                    repmat(cond_name_list(cond_idx), num_units, 1), repmat(cond_idx, num_units, 1), repmat(sess_id, num_units, 1), repmat(sess_idx, num_units, 1), (1:num_units)', cond_unit_list, ...
                    repmat(param.epoch_names(epoch_idx), num_units, 1), repmat(epoch_idx, num_units, 1), repmat(num_trials_RR, num_units, 1), repmat(num_trials_LL, num_units, 1), ...
                    prep_vals_RR(:, prep_idx), prep_vals_LL(:, prep_idx), abs(prep_vals_RR(:, prep_idx)-prep_vals_LL(:, prep_idx)), ...
                    'VariableNames', {'Preprocessing', 'PreprocessingLabel', 'ValueUnit', 'Condition', 'ConditionIndex', 'SessionID', 'Session', 'SessionUnit', ...
                    'ConditionUnit', 'Epoch', 'EpochIndex', 'NTrialsRR', 'NTrialsLL', 'RRValue', 'LLValue', 'AbsRRMinusLLValue'});
                value_table = [value_table; temp_table];
            end
        end

        cond_unit_offset = cond_unit_offset + num_units;
    end
end

%% Pooled-neuron summaries and statistics
% Kruskal-Wallis (Holm across epochs), pairwise rank-sum (Bonferroni across pairs)
summary_table = Calculate_pooled_summary(value_table, param, "AbsRRMinusLLValue", param.prep_names);
[omnibus_table, pairwise_table] = Calculate_pooled_statistics(value_table, param, "AbsRRMinusLLValue", param.prep_names);

%%% Number of units not shown in the cropped [0, 2] histograms
overflow_table = table();
for prep_idx = 1:num_prep
    for epoch_idx = 1:num_epochs
        for cond_idx = 1:numel(cond_name_list)
            selected_idx = value_table.Preprocessing == param.prep_names(prep_idx) & value_table.EpochIndex == epoch_idx & value_table.ConditionIndex == cond_idx;
            selected_vals = value_table.AbsRRMinusLLValue(selected_idx);
            selected_vals = selected_vals(isfinite(selected_vals));
            num_total = numel(selected_vals);
            num_excluded = sum(selected_vals > param.abs_diff_xlim(2));

            temp_table = table(param.prep_names(prep_idx), param.prep_labels(prep_idx), param.epoch_names(epoch_idx), epoch_idx, cond_name_list(cond_idx), cond_idx, ...
                param.abs_diff_xlim(2), num_excluded, num_total, 100*num_excluded/num_total, ...
                'VariableNames', {'Preprocessing', 'PreprocessingLabel', 'Epoch', 'EpochIndex', 'Condition', 'ConditionIndex', 'UpperXLimit', 'NExcludedAboveLimit', 'NTotal', 'PercentExcluded'});
            overflow_table = [overflow_table; temp_table];
        end
    end
end

%%% RR and LL of the no-smoothing z-score
no_smoothing_name = "TrialAvgZScore";
rr_ns_summary_table = Calculate_pooled_summary(value_table, param, "RRValue", no_smoothing_name);
[rr_ns_omnibus_table, rr_ns_pairwise_table] = Calculate_pooled_statistics(value_table, param, "RRValue", no_smoothing_name);
ll_ns_summary_table = Calculate_pooled_summary(value_table, param, "LLValue", no_smoothing_name);
[ll_ns_omnibus_table, ll_ns_pairwise_table] = Calculate_pooled_statistics(value_table, param, "LLValue", no_smoothing_name);

%%% Before vs Atropine only (one pair per epoch, so no Bonferroni correction)
before_atropine_param = param;
before_atropine_param.cond_name_list = param.cond_name_list(1:2);
before_atropine_param.cond_colors = param.cond_colors(1:2, :);
before_atropine_param.cond_pairs = [1, 2];
before_atropine_param.bar_cond_order = [1, 2];

before_atropine_values = value_table(value_table.Preprocessing == no_smoothing_name & ismember(value_table.ConditionIndex, [1, 2]), :);
before_atropine_summary_table = Calculate_pooled_summary(before_atropine_values, before_atropine_param, "AbsRRMinusLLValue", no_smoothing_name);
[before_atropine_omnibus_table, before_atropine_pairwise_table] = Calculate_pooled_statistics(before_atropine_values, before_atropine_param, "AbsRRMinusLLValue", no_smoothing_name);

%%% Print Kruskal-Wallis p-values
for prep_idx = 1:num_prep
    fprintf('\nPreprocessing: %s\n', param.prep_labels(prep_idx));
    for epoch_idx = 1:num_epochs
        selected_idx = omnibus_table.Preprocessing == param.prep_names(prep_idx) & omnibus_table.EpochIndex == epoch_idx;
        fprintf('%s: KW p = %.3g (Holm %.3g)\n', param.epoch_names(epoch_idx), omnibus_table.PRaw(selected_idx), omnibus_table.PHolmAcrossEpochs(selected_idx));
    end
end

%% Save results
fpath = strcat(param.save_path, param.version, parts_path);
if ~isfolder(fpath)
    mkdir(fpath);
end
fname = strcat(fpath, 'Fig5Di3_', param.version, '_AbsRRminusLL');

save(strcat(fname, '.mat'), 'baseline_table', 'value_table', 'summary_table', 'omnibus_table', 'pairwise_table', 'overflow_table', 'rr_ns_summary_table', 'rr_ns_omnibus_table', 'rr_ns_pairwise_table', 'll_ns_summary_table', 'll_ns_omnibus_table', 'll_ns_pairwise_table', 'before_atropine_summary_table', 'before_atropine_omnibus_table', 'before_atropine_pairwise_table', 'before_atropine_param', 'param', '-v7.3')
load(strcat(fname, '.mat'))

writetable(baseline_table, strcat(fname, '_BaselineParameters.csv'));
writetable(value_table, strcat(fname, '_NeuronEpochValues.csv'));
writetable(summary_table, strcat(fname, '_PooledSummary.csv'));
writetable(omnibus_table, strcat(fname, '_OmnibusStats.csv'));
writetable(pairwise_table, strcat(fname, '_PairwiseStats.csv'));
writetable(overflow_table, strcat(fname, '_X0to2_ExcludedCounts.csv'));
writetable(rr_ns_summary_table, strcat(fname, '_NoSmoothing_RR_PooledSummary.csv'));
writetable(rr_ns_omnibus_table, strcat(fname, '_NoSmoothing_RR_OmnibusStats.csv'));
writetable(rr_ns_pairwise_table, strcat(fname, '_NoSmoothing_RR_PairwiseStats.csv'));
writetable(ll_ns_summary_table, strcat(fname, '_NoSmoothing_LL_PooledSummary.csv'));
writetable(ll_ns_omnibus_table, strcat(fname, '_NoSmoothing_LL_OmnibusStats.csv'));
writetable(ll_ns_pairwise_table, strcat(fname, '_NoSmoothing_LL_PairwiseStats.csv'));
writetable(before_atropine_summary_table, strcat(fname, '_TrialAvgZScore_BeforeAtropine_PooledSummary.csv'));
writetable(before_atropine_omnibus_table, strcat(fname, '_TrialAvgZScore_BeforeAtropine_OmnibusStats.csv'));
writetable(before_atropine_pairwise_table, strcat(fname, '_TrialAvgZScore_BeforeAtropine_PairwiseStats.csv'));

%% Plot each row of the TrialAvgZScore [0, 2] figure (Before vs Atropine)
% Rows: mean +/- SEM bar plot, probability histogram, cumulative histogram
row_file_tags = {'Row1_Bar', 'Row2_Histogram', 'Row3_CumulativeHistogram'};
for row_idx = 1:3
    fig = Plot_distribution_figure(before_atropine_values, before_atropine_summary_table, before_atropine_pairwise_table, no_smoothing_name, "AbsRRMinusLLValue", "|RR - LL|", true, true, before_atropine_param, row_idx);
    set(fig, 'Name', ['Before vs Atropine - ', row_file_tags{row_idx}]);
    print(fig, strcat(fname, '_TrialAvgZScore_X0to2_BeforeAtropine_', row_file_tags{row_idx}, '.svg'), '-dsvg', '-painters');
end
