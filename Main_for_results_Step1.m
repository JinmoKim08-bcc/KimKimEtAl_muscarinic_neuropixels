clc
clear
close all

%% Hyperparameters to adjust
param = struct();

% add path NPY library for reading Phy
npy_path = "utils\npy-matlab\npy-matlab\";
addpath(npy_path)

%% Hyperparameters that most do not need to be adjusted
param.fs = 30000; % [Hz]; sampling rate
param.bin_time = 0.001; % [sec]; temporal binning size to calculate firing rate
param.smoothing_pts = 100; % [idx]; smoothing data points
param.boxcar.a = 1;
param.boxcar.b = ones(1, param.smoothing_pts);
param.ITI_time = 2; % [sec]; Inter-trial interval duration
param.sample_time = 0.1; % [sec]; timing given a sample
param.delay_time = 2; % [sec]
param.gocue_time = param.sample_time + param.delay_time; % [sec]
param.response_dur = 0.5; % [sec]
param.t_max = param.sample_time+ param.delay_time + param.response_dur; % [sec]; maximal end time to analyze
param.t_max_truncated = param.sample_time+ param.delay_time + param.response_dur; % [sec];
param.epochs.stimulus = [0, param.sample_time];
param.epochs.choice = [param.gocue_time-0.6, param.gocue_time];
param.epochs.action = [param.gocue_time, param.gocue_time+0.6];
param.epochs.outcome = [param.gocue_time+2, param.gocue_time+3];
param.eps = 1e-1;
param.time_epochs = [param.sample_time, (param.sample_time+param.delay_time), (param.sample_time+param.delay_time+param.response_dur)];
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

save_path = 'C:\Users\sirom\Documents\KimKimEtAl_muscarnic_neuropixels\';

param.num_iter = 100; % number of shuffling

param.version = "v16fd";

param.do_norm = 4; % 0: do not; 1: z-score; 2: normalized by baseline mean == x / mu; 3: mean subtraction; 4: z-score without param.eps; 5: mean subtraction and normalized by its norm

%% Read data and calculate firing vector
% extracted results from Main_for_Neuropixel_Kilosort_CD_v14g_for_ext_for_shuf.m
cond_name_list= ["Before"; "Atropine"; "Saline"];
% cond_name_list= ["Before"];

num_sess_list = [5, 4, 3];

max_pts = param.t_max_truncated*1/param.bin_time; %
ITI_dur = (1:param.ITI_time*param.fs)'; % 2 [sec]
ITI_dur_bin = numel(Calculate_bins_v2(ITI_dur, param.fs, param.bin_time));

elim_results = struct();

results = struct();
for cond_idx = 1:numel(cond_name_list)
    elim_results(cond_idx).elim = [];
%     elim_results(cond_idx).orig = [];

    for sess_idx = 1:num_sess_list(cond_idx)
        fprintf('cond: [%d/%d], sess: [%d/%d]\n', cond_idx, numel(cond_name_list), sess_idx, num_sess_list(cond_idx));
        load(strcat(param.save_path, 'r_', param.version, '_', cond_name_list(cond_idx), '_sess', num2str(sess_idx), '.mat'))
        
        %%% Find significant units %%%
        accum_elim_indices = [];
        for trial_idx = 1:numel(trial_types)
            selected_trial_type = trial_types(trial_idx);

            % obtain PSTH
            temp_spks_all = r.(selected_trial_type).spks / param.bin_time;
            results(cond_idx).sess(sess_idx).(selected_trial_type).allTrials.vals = temp_spks_all;

            % trial-averaged PSTH
            temp_spks = sum(temp_spks_all, 3) / size(temp_spks_all, 3);
            results(cond_idx).sess(sess_idx).(selected_trial_type).vals = temp_spks;

            % z-score normalization of firing rates by the baseline for each trial
            temp_spks_all = Normalize_firing_rates_TrialWise(temp_spks_all, param);
            temp_spks = Normalize_firing_rates(temp_spks, param);

            results(cond_idx).sess(sess_idx).(selected_trial_type).norm.vals = temp_spks;
                
            % temporal smoothing of the firing rates
            results(cond_idx).sess(sess_idx).(selected_trial_type).smooth.vals = movmean(results(cond_idx).sess(sess_idx).(selected_trial_type).norm.vals(:, ITI_dur_bin+1:(ITI_dur_bin+param.t_max*param.fs/(param.bin_time*param.fs))), param.smoothing_pts, 2);
 
            %%% Find erroneous units
            accum_elim_indices = Find_erroneous_units(results, cond_idx, sess_idx, selected_trial_type, accum_elim_indices, param);
        end

        accum_elim_indices = unique(accum_elim_indices);
        fprintf("Cond: %d, Sess: %d\n", cond_idx, sess_idx);
        results = Eliminate_erroneous_units_v2(results, trial_types, cond_idx, sess_idx, accum_elim_indices);

        % PSTH for correct trials
        temp_x = (results(cond_idx).sess(sess_idx).RR.vals + results(cond_idx).sess(sess_idx).LL.vals) / 2; % PSTH for correct trials
        temp_x = Normalize_firing_rates(temp_x, param);
        results(cond_idx).sess(sess_idx).cor.norm.vals = temp_x;
        x_baseline = movmean(results(cond_idx).sess(sess_idx).cor.norm.vals(:, 1:ITI_dur_bin), param.smoothing_pts, 2);
        x = movmean(results(cond_idx).sess(sess_idx).cor.norm.vals(:, ITI_dur_bin+1:(ITI_dur_bin+max_pts)), param.smoothing_pts, 2);
        results(cond_idx).sess(sess_idx).cor.smooth.vals = x;

          time_bin = (0:param.bin_time:(param.t_max-param.bin_time))';
        figure
        hold on
        plot(time_bin, mean(x, 1), 'k')
%         plot(time_bin, mean(x_shuf, 1), 'r')
        xlabel('Time [sec]')
        ylabel('Normalized firing rate')
        title('Normalized and smoothed')
        set(gca, 'FontSize', 15)

        if numel(x_baseline) ~= 0
            % Find significantly responsive neurons
            [sig_idx, latency] = Find_response_latency(x_baseline, x, max_pts);
            results(cond_idx).sess(sess_idx).sig_idx = sig_idx;
            results(cond_idx).sess(sess_idx).latency = latency;
        else
            results(cond_idx).sess(sess_idx).sig_idx = [];
            results(cond_idx).sess(sess_idx).latency = [];
        end
        
        %%% Rasters
        temp_unit = r.RR.unit;
        temp_unit(accum_elim_indices) = [];
        results(cond_idx).sess(sess_idx).RR.unit = temp_unit;
        temp_unit = r.RL.unit;
        temp_unit(accum_elim_indices) = [];
        results(cond_idx).sess(sess_idx).RL.unit = temp_unit;
        temp_unit = r.LL.unit;
        temp_unit(accum_elim_indices) = [];
        results(cond_idx).sess(sess_idx).LL.unit = temp_unit;
        temp_unit = r.LR.unit;
        temp_unit(accum_elim_indices) = [];
        results(cond_idx).sess(sess_idx).LR.unit = temp_unit;

        elim_results(cond_idx).elim = [elim_results(cond_idx).elim; accum_elim_indices];
%         elim_results(cond_idx).orig = [elim_results(cond_idx).orig; ];

        results(cond_idx).sess(sess_idx).RR.trial_info = r.RR.trial_info;
        results(cond_idx).sess(sess_idx).RL.trial_info = r.RL.trial_info;
        results(cond_idx).sess(sess_idx).LL.trial_info = r.LL.trial_info;
        results(cond_idx).sess(sess_idx).LR.trial_info = r.LR.trial_info;
        
    end
end
save(strcat(param.save_path, '\results_', param.version, '.mat'), 'results', '-v7.3')
