clc
clear
close all

%% Hyperparameters to adjust
param = struct();

% add path NPY library for reading Phy
npy_path = "utils\npy-matlab\npy-matlab\";
addpath(npy_path)

save_path = 'C:\Users\sirom\Documents\KimKimEtAl_muscarnic_neuropixels\';
parts_path = '\fig4\';

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
param.num_bootstrap = 20;

%% Read calculated firing vectors
cond_name_list= ["Before"; "Atropine"; "Saline"];

num_sess_list = [5, 4, 3];

max_pts = param.t_max_truncated*1/param.bin_time; %
ITI_dur = (1:param.ITI_time*param.fs)'; % 2 [sec]
ITI_dur_bin = numel(Calculate_bins_v2(ITI_dur, param.fs, param.bin_time));

load(strcat(param.save_path, 'results_', param.version, '.mat'));

%% Heatmap plot
clc
close all

fpath = strcat(param.save_path, '\');

sep_trials = ["cor"; "CR"; "CL"; "CR-CL"; "incor"; "ER"; "EL"; "ER-EL"; "CR-ER"; "CL-EL"; "Cor+Incor"];

time_bin = -param.ITI_time:param.bin_time:param.t_max_truncated-param.bin_time;
% time_bin = (0:param.bin_time:param.t_max-param.bin_time)';
z_score_bounds = [0, 3];
time_bounds = [-1, time_bin(end)];

stats = struct();

for cond_idx = 1 % :numel(cond_name_list)
    CR = [];
    CL = [];
    ER = [];
    EL = [];

    for sess_idx = 1:num_sess_list(cond_idx)
        CR = [CR; results(cond_idx).sess(sess_idx).RR.norm.vals];
        CL = [CL; results(cond_idx).sess(sess_idx).LL.norm.vals];
        ER = [ER; results(cond_idx).sess(sess_idx).LR.norm.vals];
        EL = [EL; results(cond_idx).sess(sess_idx).RL.norm.vals];
    end

    for trial_idx = [1, 5] % [1, 2, 3, 5, 6, 7] %1:numel(sep_trials)
        selected_trial_type = sep_trials(trial_idx);

        x_baseline = [];
        x = [];

        switch trial_idx
            case 1
                temp_spks = (CR+CL)/2;
            case 2
                temp_spks = CR;
            case 3
                temp_spks = CL;
            case 4
                temp_spks = (CR-CL)/2;
            case 5
                temp_spks = (ER+EL)/2;
            case 6
                temp_spks = ER;
            case 7
                temp_spks = EL;
            case 8
                temp_spks = (ER-EL)/2;
            case 9
                temp_spks = (CR-ER)/2;
            case 10
                temp_spks = (CL-EL)/2;
            case 11
                temp_spks = (CR+CL+ER+EL)/4;
        end

        x_baseline = movmean(temp_spks(:, 1:ITI_dur_bin), param.smoothing_pts, 2);
        x = movmean(temp_spks(:, (ITI_dur_bin+1):end), param.smoothing_pts, 2);

        if trial_idx == 1
            [sig_idx, latency] = Find_response_latency(x_baseline, x, max_pts);

            [~, indices] = sort(latency, 'ascend');

            stats(cond_idx).sig_idx = sig_idx;
            stats(cond_idx).latency = latency;
        end

        x = movmean(temp_spks, param.smoothing_pts, 2);

        y = 1:size(sig_idx, 1); % unit indices

        fig = figure;
        h = imagesc(time_bin, y, x(sig_idx(indices), :));
        hold on
        xline(0, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
        xline(param.sample_time, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
        xline(param.gocue_time, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
        colormap turbo
        colorbar
        caxis(z_score_bounds)

        xlabel('Time [sec]')
        ylabel('Units')
        xlim(time_bounds)
        set(gca,'FontSize',15)
        title(strcat(cond_name_list(cond_idx), '-', selected_trial_type))
        print(fig, strcat(fpath, parts_path, '\', cond_name_list(cond_idx), '-', selected_trial_type), '-dsvg')
    end
end

% %% Firing rates plot
% for cond_idx = 1 % :numel(cond_name_list)
%     CR = [];
%     CL = [];
%     ER = [];
%     EL = [];
% 
%     for sess_idx = 1:num_sess_list(cond_idx)
%         CR = [CR; results(cond_idx).sess(sess_idx).RR.norm.vals];
%         CL = [CL; results(cond_idx).sess(sess_idx).LL.norm.vals];
%         ER = [ER; results(cond_idx).sess(sess_idx).LR.norm.vals];
%         EL = [EL; results(cond_idx).sess(sess_idx).RL.norm.vals];
%     end
% 
%     for trial_idx = [1] % , 5] % [1, 2, 3, 5, 6, 7] %1:numel(sep_trials)
%         selected_trial_type = sep_trials(trial_idx);
% 
%         x_baseline = [];
%         x = [];
% 
%         switch trial_idx
%             case 1
%                 temp_spks = (CR+CL)/2;
%             case 2
%                 temp_spks = CR;
%             case 3
%                 temp_spks = CL;
%             case 4
%                 temp_spks = (CR-CL)/2;
%             case 5
%                 temp_spks = (ER+EL)/2;
%             case 6
%                 temp_spks = ER;
%             case 7
%                 temp_spks = EL;
%             case 8
%                 temp_spks = (ER-EL)/2;
%             case 9
%                 temp_spks = (CR-ER)/2;
%             case 10
%                 temp_spks = (CL-EL)/2;
%             case 11
%                 temp_spks = (CR+CL+ER+EL)/4;
%         end
% 
%         x_baseline = movmean(temp_spks(:, 1:ITI_dur_bin), param.smoothing_pts, 2);
%         x = movmean(temp_spks(:, (ITI_dur_bin+1):end), param.smoothing_pts, 2);
% 
%         if trial_idx == 1
%             [sig_idx, latency] = Find_response_latency(x_baseline, x, max_pts);
% 
%             [~, indices] = sort(latency, 'ascend');
% 
%             stats(cond_idx).sig_idx = sig_idx;
%             stats(cond_idx).latency = latency;
%         end
% 
%         x = movmean(temp_spks, param.smoothing_pts, 2);
% 
%         y = 1:size(sig_idx, 1); % unit indices
% 
%         mycolormap = jet(numel(y));
% 
%         fig = figure;
%         % h = imagesc(time_bin, y, x(sig_idx(indices), :));
%         plot(time_bin, x')
%         hold on
%         xline(0, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
%         xline(param.sample_time, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
%         xline(param.gocue_time, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
%         colormap(mycolormap)
%         colorbar
%         % caxis(z_score_bounds)
% 
%         xlabel('Time [sec]')
%         ylabel('Firing rates')
%         xlim(time_bounds)
%         set(gca,'FontSize',15)
%         title(strcat(cond_name_list(cond_idx), '-', selected_trial_type))
%         % print(fig, strcat(fpath, cond_name_list(cond_idx), '-', selected_trial_type), '-dsvg')
%     end
% end

% %% Firing rates histogram
% for cond_idx = 1 % :numel(cond_name_list)
%     CR = [];
%     CL = [];
%     ER = [];
%     EL = [];
% 
%     for sess_idx = 1:num_sess_list(cond_idx)
%         CR = [CR; results(cond_idx).sess(sess_idx).RR.norm.vals];
%         CL = [CL; results(cond_idx).sess(sess_idx).LL.norm.vals];
%         ER = [ER; results(cond_idx).sess(sess_idx).LR.norm.vals];
%         EL = [EL; results(cond_idx).sess(sess_idx).RL.norm.vals];
%     end
% 
%     for trial_idx = [1] % , 5] % [1, 2, 3, 5, 6, 7] %1:numel(sep_trials)
%         selected_trial_type = sep_trials(trial_idx);
% 
%         x_baseline = [];
%         x = [];
% 
%         switch trial_idx
%             case 1
%                 temp_spks = (CR+CL)/2;
%             case 2
%                 temp_spks = CR;
%             case 3
%                 temp_spks = CL;
%             case 4
%                 temp_spks = (CR-CL)/2;
%             case 5
%                 temp_spks = (ER+EL)/2;
%             case 6
%                 temp_spks = ER;
%             case 7
%                 temp_spks = EL;
%             case 8
%                 temp_spks = (ER-EL)/2;
%             case 9
%                 temp_spks = (CR-ER)/2;
%             case 10
%                 temp_spks = (CL-EL)/2;
%             case 11
%                 temp_spks = (CR+CL+ER+EL)/4;
%         end
% 
%         x_baseline = movmean(temp_spks(:, 1:ITI_dur_bin), param.smoothing_pts, 2);
%         x = movmean(temp_spks(:, (ITI_dur_bin+1):end), param.smoothing_pts, 2);
% 
%         if trial_idx == 1
%             [sig_idx, latency] = Find_response_latency(x_baseline, x, max_pts);
% 
%             %             [sig_idx, latency, segments] = Find_response_latency_consecutive(...
%             %                 x_baseline, x, max_pts);
% 
%             [~, indices] = sort(latency, 'ascend');
% 
%             stats(cond_idx).sig_idx = sig_idx;
%             stats(cond_idx).latency = latency;
%             %             stats(cond_idx).segments = segments;
%         end
% 
%         x = movmean(temp_spks, param.smoothing_pts, 2);
% 
%         y = 1:size(sig_idx, 1); % unit indices
% 
%         t_lb_idx = find(time_bin == 0);
%         [~, t_ub_idx] = min(abs(time_bin - 0.1));
%         t_ub_idx = t_ub_idx - 1;
% 
%         x_stim = mean(x(:, t_lb_idx:t_ub_idx), 2);
% 
%         fig = figure;
% 
%         edges = min(x_stim):0.1:max(x_stim);
%         hist(x_stim, edges)
%         xlabel('Firing rates')
%         ylabel('#(Units)')
%         set(gca, 'FontSize', 15)
%     end
% end


%% Epoch-selectivity
param.epochs.stimulus = [0, param.sample_time];
param.epochs.delay = [param.sample_time, param.sample_time+param.delay_time];
param.epochs.response = [param.gocue_time, param.gocue_time+0.5];


epoch_names = ["stimulus"; "delay"; "response"];
selectivity_eqs = [1, 3, 2, 4;
    1, 4, 2, 3;
    1, 4, 2, 3];

offset = param.ITI_time*param.fs_bin;

epoch_time_idx = [1, param.time_epochs(1)*param.fs_bin;
                  1+param.time_epochs(1)*param.fs_bin, param.time_epochs(2)*param.fs_bin;
                  1+param.time_epochs(2)*param.fs_bin, param.time_epochs(3)*param.fs_bin];

for cond_idx = 1:numel(cond_name_list)
    num_sess = size(results(cond_idx).sess, 2);

    for sess_idx = 1:num_sess
        num_unit = size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1);

        for epoch_idx = 1:numel(epoch_names)
            fprintf('Cond_idx: %d, Sess_idx: %d, Epoch_idx: %d\n', cond_idx, sess_idx, epoch_idx)

            seqs = [];
            epoch_window = param.epochs.(epoch_names(epoch_idx));
            lb_idx = offset + floor(epoch_window(1)*param.fs_bin)+1;
            ub_idx = offset + floor(epoch_window(2)*param.fs_bin);
            %             results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).group1 = zeros(num_unit, 1);
            %             results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).group2 = zeros(num_unit, 1);

            temp_a = results(cond_idx).sess(sess_idx).(trial_types(selectivity_eqs(epoch_idx, 1))).allTrials.vals;
            temp_a(isnan(temp_a)) = 0;
            temp_a = (temp_a - mean(temp_a(:, 1:offset, :), 2)); % ./ (std(temp_a(:, 1:offset, :), [], 2) + 1e-3);
            temp_a = temp_a(:, ITI_dur_bin+1:end, :);
            temp_a = movmean(temp_a, param.smoothing_pts, 2);

            temp_b = results(cond_idx).sess(sess_idx).(trial_types(selectivity_eqs(epoch_idx, 2))).allTrials.vals;
            temp_b(isnan(temp_b)) = 0;
            temp_b = (temp_b - mean(temp_b(:, 1:offset, :), 2)); % ./ (std(temp_b(:, 1:offset, :), [], 2) + 1e-3);
            temp_b = temp_b(:, ITI_dur_bin+1:end, :);
            temp_b = movmean(temp_b, param.smoothing_pts, 2);

            temp_c = results(cond_idx).sess(sess_idx).(trial_types(selectivity_eqs(epoch_idx, 3))).allTrials.vals;
            temp_c(isnan(temp_c)) = 0;
            temp_c = (temp_c - mean(temp_c(:, 1:offset, :), 2)); % ./ (std(temp_c(:, 1:offset, :), [], 2) + 1e-3);
            temp_c = temp_c(:, ITI_dur_bin+1:end, :);
            temp_c = movmean(temp_c, param.smoothing_pts, 2);

            temp_d = results(cond_idx).sess(sess_idx).(trial_types(selectivity_eqs(epoch_idx, 4))).allTrials.vals;
            temp_d(isnan(temp_d)) = 0;
            temp_d = (temp_d - mean(temp_d(:, 1:offset, :), 2)); % ./ (std(temp_d(:, 1:offset, :), [], 2) + 1e-3);
            temp_d = temp_d(:, ITI_dur_bin+1:end, :);
            temp_d = movmean(temp_d, param.smoothing_pts, 2);

            a = squeeze(mean(temp_a(:,epoch_time_idx(epoch_idx,1):epoch_time_idx(epoch_idx,2), :), 2));
            b = squeeze(mean(temp_b(:,epoch_time_idx(epoch_idx,1):epoch_time_idx(epoch_idx,2), :), 2));
            c = squeeze(mean(temp_c(:,epoch_time_idx(epoch_idx,1):epoch_time_idx(epoch_idx,2), :), 2));
            d = squeeze(mean(temp_d(:,epoch_time_idx(epoch_idx,1):epoch_time_idx(epoch_idx,2), :), 2));

            seqs(1).val = a;
            seqs(2).val = b;
            seqs(3).val = c;
            seqs(4).val = d;

            results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).groups(1).val = [seqs(selectivity_eqs(epoch_idx, 1)).val, seqs(selectivity_eqs(epoch_idx, 2)).val];
            results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).groups(2).val = [seqs(selectivity_eqs(epoch_idx, 3)).val, seqs(selectivity_eqs(epoch_idx, 4)).val];
        end
    end
end

% selectives_neurons = struct('locs', struct());

for cond_idx = 1:numel(cond_name_list)
    num_sess = size(results(cond_idx).sess, 2);

    for sess_idx = 1:num_sess
        %         selectives_neurons(cond_idx).sess(sess_idx).stimulus = struct('a', [], 'b', [], 'c', [], 'd', []);
        %         selectives_neurons(cond_idx).sess(sess_idx).delay = struct('a', [], 'b', [], 'c', [], 'd', []);
        %         selectives_neurons(cond_idx).sess(sess_idx).response = struct('a', [], 'b', [], 'c', [], 'd', []);
        %         selectives_neurons(cond_idx).sess(sess_idx).choice = struct('a', [], 'b', [], 'c', [], 'd', []);
        %         selectives_neurons(cond_idx).sess(sess_idx).action = struct('a', [], 'b', [], 'c', [], 'd', []);
        %         selectives_neurons(cond_idx).sess(sess_idx).outcome = struct('a', [], 'b', [], 'c', [], 'd', []);

        num_unit = size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1);

        for epoch_idx = 1:numel(epoch_names)
            results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).RR = [];
            results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).LL = [];

            x = results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).groups(1).val;
            y = results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).groups(2).val;

            num_selective = 0;
            selective_unit_indices = [];
            for unit_idx = 1:num_unit
                fprintf('Cond: %d, Sess: %d, Epoch: %d, unit: %d\n', cond_idx, sess_idx, epoch_idx, unit_idx)

                [p, h] = ranksum(x(unit_idx, :), y(unit_idx, :)); %, 'alpha', 0.05);
                if h == 1
                    RR = results(cond_idx).sess(sess_idx).RR.allTrials.vals(unit_idx, :, :);
                    RR = squeeze(mean(RR, 3));
                    RR = movmean(RR, param.smoothing_pts, 2);
                    RR = (RR - mean(RR(:, 1:ITI_dur_bin), 2)) ./ (std(RR(:, 1:ITI_dur_bin), [], 2)+1e-3);

                    LL = results(cond_idx).sess(sess_idx).LL.allTrials.vals(unit_idx, :, :);
                    LL = squeeze(mean(LL, 3));
                    LL = movmean(LL, param.smoothing_pts, 2);
                    LL = (LL - mean(LL(:, 1:ITI_dur_bin), 2)) ./ (std(LL(:, 1:ITI_dur_bin), [], 2)+1e-3);

                    results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).RR = [results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).RR; RR];
                    results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).LL = [results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).LL; LL];

                    num_selective = num_selective + 1;
                    selective_unit_indices = [selective_unit_indices; unit_idx];
                end
            end
            results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).num_unit = num_unit;
            results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).num_selective = num_selective;
            results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).selective_unit_indices = selective_unit_indices;
        end
    end
end

% %%% number of selectivity
% for cond_idx = 1:numel(cond_name_list)
%     results(cond_idx).traces = struct();
%     for epoch_idx = 1:numel(epoch_names)
%         num_total_unit = 0;
%         num_total_selective_units = 0;
%
%         results(cond_idx).traces.(epoch_names(epoch_idx)).RR = [];
%         results(cond_idx).traces.(epoch_names(epoch_idx)).LL = [];
%
%         sess_wise_units = [];
%
%         for sess_idx = 1:num_sess
%             num_total_unit = num_total_unit + results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).num_unit;
%             num_total_selective_units = num_total_selective_units + results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).num_selective;
%
%             selected_indices = results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).selective_unit_indices;
%
%             temp_RR_traces = results(cond_idx).sess(sess_idx).RR.norm.vals(selected_indices, :);
%             temp_RR_traces = movmean(temp_RR_traces, param.smoothing_pts, 2);
%             temp_LL_traces = results(cond_idx).sess(sess_idx).LL.norm.vals(selected_indices, :);
%             temp_LL_traces = movmean(temp_LL_traces, param.smoothing_pts, 2);
%
%             results(cond_idx).traces.(epoch_names(epoch_idx)).RR = [results(cond_idx).traces.(epoch_names(epoch_idx)).RR; temp_RR_traces];
%             results(cond_idx).traces.(epoch_names(epoch_idx)).LL = [results(cond_idx).traces.(epoch_names(epoch_idx)).LL; temp_LL_traces];
%
%             sess_wise_units = [sess_wise_units; num_total_selective_units];
%         end
%     end
% end

%% Epoch-selectivity: firing plot for all selective units
% clc
% close all
% 
% for cond_idx = 1:numel(cond_name_list)
%     num_sess = size(results(cond_idx).sess, 2);
% 
%     for sess_idx = 1:num_sess
%         for epoch_idx = 1:numel(epoch_names)
%             num_unit = size(results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).RR, 1);
%             time_bin = -param.ITI_time:param.bin_time:param.t_max_truncated-param.bin_time;
% 
%             for unit_idx = 1:num_unit
%                 fig = figure;
%                 hold on
%                 plot(time_bin, results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).RR(unit_idx, :), 'LineWidth', 1.5, 'color', 'b')
%                 plot(time_bin, results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).LL(unit_idx, :), 'LineWidth', 1.5, 'color', 'r')
%                 xline(0, ':', 'Color', [0, 0, 0], 'LineWidth', 1.5)
%                 xline(param.sample_time, ':', 'Color', [0, 0, 0], 'LineWidth', 1.5)
%                 xline(param.gocue_time, ':', 'Color', [0, 0, 0], 'LineWidth', 1.5)
%                 xlim(time_bounds)
%                 xlabel('Time [sec]')
%                 ylabel('Normalized firing rate [a.u.]')
%                 set(gca, 'FontSize', 15)
%                 print(fig, strcat(fpath, parts_path, '\', '\traces\', epoch_names(epoch_idx), '\cond', num2str(cond_idx), '-sess', num2str(sess_idx), '-unit', num2str(unit_idx)), '-dpng')
%                 close all
%             end
%         end
%     end
% end

%% Epoch-selectivity: raster plot for representative units
%%% stimulus selective unit
sess_idx = 3;
cond_idx = 3;
epoch_idx = 1;
% sess_idx = 1;
% cond_idx = 3;
% epoch_idx = 3;
selective_indices = results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).selective_unit_indices;
selected_unit_idx = selective_indices(4);
%%%% RR
temp_RR = results(cond_idx).sess(sess_idx).RR.unit(selected_unit_idx).rasters;
%%%% LL
temp_LL = results(cond_idx).sess(sess_idx).LL.unit(selected_unit_idx).rasters;

% temp_RR, temp_LL: N×K cell array
[~, K_RR] = size(temp_RR);
[~, K_LL] = size(temp_LL);

fig = figure('Units', 'pixels', 'Position', [753, 142, 577, 748]);
set(fig, 'Renderer', 'painters', 'GraphicsSmoothing', 'off');

subplot(3, 1, [1, 2])
hold on;

time_bin = (-param.ITI_time:param.bin_time:param.t_max-param.bin_time)';

% ——— RR (blue) ———
for tr = 1:K_RR
    % Concatenate all row data of the same trial (tr) into one vector
    tRR = vertcat( temp_RR{:, tr} );
    yRR = tr + K_LL;                    % offset into the upper half
    scatter( tRR, yRR*ones(size(tRR)), 12, 'b', 'filled' );
end

% ——— LL (red) ———
for tr = 1:K_LL
    tLL = vertcat( temp_LL{:, tr} );
    yLL = tr;                        % lower half
    scatter( tLL, yLL*ones(size(tLL)), 12, 'r', 'filled' );
end

xline(0, ':', 'Color', 'k', 'LineWidth', 3)
xline(param.sample_time, ':', 'Color', 'k', 'LineWidth', 3)
xline(param.gocue_time, ':', 'Color', 'k', 'LineWidth', 3)

% Axes and labels
xlabel('Time [sec]');
ylabel('Trials');
ylim([1, K_RR+K_LL])
xlim(time_bounds)
set(gca, 'FontSize', 15)

b = ones(1,param.smoothing_pts)/param.smoothing_pts;

subplot(3, 1, 3)
hold on;
plot(time_bin, movmean(results(cond_idx).sess(sess_idx).RR.norm.vals(selected_unit_idx, :), param.smoothing_pts), 'LineWidth', 1.5, 'color', 'b')
plot(time_bin, movmean(results(cond_idx).sess(sess_idx).LL.norm.vals(selected_unit_idx, :), param.smoothing_pts), 'LineWidth', 1.5, 'color', 'r')
xline(0, ':', 'Color', 'k', 'LineWidth', 3)
xline(param.sample_time, ':', 'Color', 'k', 'LineWidth', 3)
xline(param.gocue_time, ':', 'Color', 'k', 'LineWidth', 3)
xlabel('Time [sec]');
ylabel('Normalized firing rate [a.u.]')
xlim(time_bounds)
set(gca, 'FontSize', 15)
print(fig, strcat(fpath, parts_path, '\', '\Raster-Stim_selective'), '-dsvg', '-painters')

%%% Delay selective unit
sess_idx = 1;
cond_idx = 1;
epoch_idx = 2;
selective_indices = results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).selective_unit_indices;
selected_unit_idx = selective_indices(46);
%%%% RR
temp_RR = results(1).sess(sess_idx).RR.unit(selected_unit_idx).rasters;
%%%% LL
temp_LL = results(1).sess(sess_idx).LL.unit(selected_unit_idx).rasters;

% temp_RR, temp_LL: N×K cell array
[~, K_RR] = size(temp_RR);
[~, K_LL] = size(temp_LL);

fig = figure('Units', 'pixels', 'Position', [753, 142, 577, 748]);
set(fig, 'Renderer', 'painters', 'GraphicsSmoothing', 'off');
subplot(3, 1, [1, 2])
hold on;

time_bin = (-param.ITI_time:param.bin_time:param.t_max-param.bin_time)';

% ——— RR (blue) ———
for tr = 1:K_RR
    % Concatenate all row data of the same trial (tr) into one vector
    tRR = vertcat( temp_RR{:, tr} );
    yRR = tr + K_LL;                    % offset into the upper half
    scatter( tRR, yRR*ones(size(tRR)), 12, 'b', 'filled' );
end

% ——— LL (red) ———
for tr = 1:K_LL
    tLL = vertcat( temp_LL{:, tr} );
    yLL = tr;                        % lower half
    scatter( tLL, yLL*ones(size(tLL)), 12, 'r', 'filled' );
end

xline(0, ':', 'Color', 'k', 'LineWidth', 3)
xline(param.sample_time, ':', 'Color', 'k', 'LineWidth', 3)
xline(param.gocue_time, ':', 'Color', 'k', 'LineWidth', 3)

% Axes and labels
xlabel('Time [sec]');
ylabel('Trials');
ylim([1, K_RR+K_LL])
xlim(time_bounds)
set(gca, 'FontSize', 15)

subplot(3, 1, 3)
hold on;
plot(time_bin, movmean(results(1).sess(sess_idx).RR.norm.vals(selected_unit_idx, :), param.smoothing_pts), 'LineWidth', 1.5, 'color', 'b')
plot(time_bin, movmean(results(1).sess(sess_idx).LL.norm.vals(selected_unit_idx, :), param.smoothing_pts), 'LineWidth', 1.5, 'color', 'r')
xline(0, ':', 'Color', 'k', 'LineWidth', 3)
xline(param.sample_time, ':', 'Color', 'k', 'LineWidth', 3)
xline(param.gocue_time, ':', 'Color', 'k', 'LineWidth', 3)
xlabel('Time [sec]');
ylabel('Normalized firing rate [a.u.]')
xlim(time_bounds)
set(gca, 'FontSize', 15)
print(fig, strcat(fpath, parts_path, '\', '\Raster-Delay_selective'), '-dsvg', '-painters')

%%% Response selective unit
% clc
% close all
sess_idx = 1;
cond_idx = 1;
epoch_idx = 3;
selective_indices = results(cond_idx).sess(sess_idx).(epoch_names(epoch_idx)).selective_unit_indices;
selected_unit_idx = selective_indices(20);
%%%% RR
temp_RR = results(1).sess(sess_idx).RR.unit(selected_unit_idx).rasters;
%%%% LL
temp_LL = results(1).sess(sess_idx).LL.unit(selected_unit_idx).rasters;

% temp_RR, temp_LL: N×K cell array
[~, K_RR] = size(temp_RR);
[~, K_LL] = size(temp_LL);

fig = figure('Units', 'pixels', 'Position', [753, 142, 577, 748]);
set(fig, 'Renderer', 'painters', 'GraphicsSmoothing', 'off');
subplot(3, 1, [1, 2])
hold on;

time_bin = (-param.ITI_time:param.bin_time:param.t_max-param.bin_time)';

% ——— RR (blue) ———
for tr = 1:K_RR
    % Concatenate all row data of the same trial (tr) into one vector
    tRR = vertcat( temp_RR{:, tr} );
    yRR = tr + K_LL;                    % offset into the upper half
    scatter( tRR, yRR*ones(size(tRR)), 12, 'b', 'filled' );
end

% ——— LL (red) ———
for tr = 1:K_LL
    tLL = vertcat( temp_LL{:, tr} );
    yLL = tr;                        % lower half
    scatter( tLL, yLL*ones(size(tLL)), 12, 'r', 'filled' );
end

xline(0, ':', 'Color', 'k', 'LineWidth', 3)
xline(param.sample_time, ':', 'Color', 'k', 'LineWidth', 3)
xline(param.gocue_time, ':', 'Color', 'k', 'LineWidth', 3)

% Axes and labels
xlabel('Time [sec]');
ylabel('Trials');
ylim([1, K_RR+K_LL])
xlim(time_bounds)
set(gca, 'FontSize', 15)

subplot(3, 1, 3)
hold on;
plot(time_bin, movmean(results(1).sess(sess_idx).RR.norm.vals(selected_unit_idx, :), param.smoothing_pts), 'LineWidth', 1.5, 'color', 'b')
plot(time_bin, movmean(results(1).sess(sess_idx).LL.norm.vals(selected_unit_idx, :), param.smoothing_pts), 'LineWidth', 1.5, 'color', 'r')
xline(0, ':', 'Color', 'k', 'LineWidth', 3)
xline(param.sample_time, ':', 'Color', 'k', 'LineWidth', 3)
xline(param.gocue_time, ':', 'Color', 'k', 'LineWidth', 3)
xlabel('Time [sec]');
ylabel('Normalized firing rate [a.u.]')
xlim(time_bounds)
set(gca, 'FontSize', 15)
print(fig, strcat(fpath, parts_path, '\', '\Raster-Response_selective'), '-dsvg', '-painters')

%% LDA analysis w/o/ time resolution for each unit (trial-averaged firing rates w/o/ LDA)
clc
close all

time_bin = -param.ITI_time:param.bin_time:param.t_max_truncated-param.bin_time;

% parameters for pseudo population
param.num_repeat = 100;
param.ratio_train = 0.8;
param.ratio_test = 1 - param.ratio_train;

selectivity_index = struct();

for cond_idx = 1:numel(cond_name_list)

    datasets = struct();
    trials = struct();
    cv_opts = struct();
    for sess_idx = 1:num_sess_list(cond_idx)
        for trial_idx = 1:numel(trial_types)
            selected_trial_type = trial_types(trial_idx);

            temp_spks = results(cond_idx).sess(sess_idx).(selected_trial_type).allTrials.vals(:, 1:numel(time_bin), :);
            num_trials = size(temp_spks, 3);
            for repeat_idx = 1:param.num_repeat
                rng(repeat_idx)

                % all_indices = 1:num_trials;
                % train_indices = randi(num_trials, 1, floor(num_trials*param.ratio_train));
                % test_indices = setdiff(all_indices, train_indices);
                perm = randperm(num_trials);
                ntr = floor(num_trials * param.ratio_train);
                train_indices = perm(1:ntr);
                test_indices  = perm(ntr+1:end);

                trials(sess_idx).(selected_trial_type).repeat(repeat_idx).testSets = test_indices;
                trials(sess_idx).(selected_trial_type).repeat(repeat_idx).trainSets = train_indices;

            end
        end
    end

    % Find the number of units
    num_units = 0;
    num_units_list = [];
    for sess_idx = 1:num_sess_list(cond_idx)
        num_units = num_units + size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1);
        num_units_list = [num_units_list; size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1)];
    end
    num_units_cumsum2 = cumsum(num_units_list);
    num_units_cumsum1 = [1; num_units_cumsum2(1:end-1)+1];
    num_units_cumsum = [num_units_cumsum1, num_units_cumsum2];

    %
    tic

    remained_unit_idx = 0;
    abnormality_check = 1;

    for unit_idx = 1:num_units

        sess_idx = find(unit_idx >= num_units_cumsum(:,1) & unit_idx <= num_units_cumsum(:,2), 1, 'first');
        local_unit_idx = unit_idx - num_units_cumsum(sess_idx,1) + 1;

        % Load RR / LL trial matrices (T_full x nTrials)
        RR = squeeze(results(cond_idx).sess(sess_idx).RR.allTrials.vals(local_unit_idx, :, :)); % [T_full x nR]
        LL = squeeze(results(cond_idx).sess(sess_idx).LL.allTrials.vals(local_unit_idx, :, :)); % [T_full x nL]
        RL = squeeze(results(cond_idx).sess(sess_idx).RL.allTrials.vals(local_unit_idx, :, :)); % [T_full x nR]
        LR = squeeze(results(cond_idx).sess(sess_idx).LR.allTrials.vals(local_unit_idx, :, :)); % [T_full x nR]

        RR = Perform_smoothing_fr(RR', param); % (nTrials x T_full)
        LL = Perform_smoothing_fr(LL', param);
        RL = Perform_smoothing_fr(RL', param); % (nTrials x T_full)
        LR = Perform_smoothing_fr(LR', param);

        RR = (RR - mean(RR(:, 1:ITI_dur_bin), 2));
        LL = (LL - mean(LL(:, 1:ITI_dur_bin), 2));
        RL = (RL - mean(RL(:, 1:ITI_dur_bin), 2));
        LR = (LR - mean(LR(:, 1:ITI_dur_bin), 2));

        proj_cor_R = [];
        proj_cor_L = [];
        proj_incor_R = [];
        proj_incor_L = [];

        acc_cor = [];
        auroc_cor = [];
        acc_incor = [];
        auroc_incor = [];

        fpr_list_cor = {};
        tpr_list_cor = {};

        fpr_list_incor = {};
        tpr_list_incor = {};

        for repeat_idx = 1:param.num_repeat
            fprintf("Cond: [%d/%d] | Unit: [%d/%d], Repeat: [%d/%d]\n", cond_idx, numel(cond_name_list), unit_idx, num_units, repeat_idx, param.num_repeat);

            %%% Correct trials
            selected_trial_type = "RR";
            train_indices = trials(sess_idx).(selected_trial_type).repeat(repeat_idx).trainSets;
            test_indices = trials(sess_idx).(selected_trial_type).repeat(repeat_idx).testSets;
            X_train_R = mean(RR(train_indices, 2001:end), 1); % 1 x T
            X_test_R = mean(RR(test_indices, 2001:end), 1); % 1 x T

            selected_trial_type = "LL";
            train_indices = trials(sess_idx).(selected_trial_type).repeat(repeat_idx).trainSets;
            test_indices = trials(sess_idx).(selected_trial_type).repeat(repeat_idx).testSets;
            X_train_L = mean(LL(train_indices, 2001:end), 1); % 1 x T
            X_test_L = mean(LL(test_indices, 2001:end), 1); % 1 x T

            [N1, T1] = size(X_train_R);
            [N2, T2] = size(X_train_L);

            X_train = [X_train_R'; X_train_L']; % (T1+T2) x 1
            Y_train = [repmat("Right", T1, 1) ; ...
                       repmat("Left" , T2, 1)]; % 2T x 1

            %%% metrics w/o/ LDA
            if mean(X_train_R) > mean(X_train_L)
                posLabel = "Right";
            else
                posLabel = "Left";
            end

            scores = X_train;
            [fpr,tpr,~,auROC] = perfcurve(Y_train, X_train, posLabel);
            auroc_cor = [auroc_cor; auROC];

            fpr_list_cor{repeat_idx} = fpr;
            tpr_list_cor{repeat_idx} = tpr;

            %%% Incorrect trials
            % selected_trial_type = "LR";
            selected_trial_type = "LR";
            train_indices = trials(sess_idx).(selected_trial_type).repeat(repeat_idx).trainSets;
            test_indices = trials(sess_idx).(selected_trial_type).repeat(repeat_idx).testSets;
            X_train_R = mean(LR(train_indices, 2001:end), 1); % 1 x T
            X_test_R = mean(LR(test_indices, 2001:end), 1); % 1 x T

            selected_trial_type = "RL";
            train_indices = trials(sess_idx).(selected_trial_type).repeat(repeat_idx).trainSets;
            test_indices = trials(sess_idx).(selected_trial_type).repeat(repeat_idx).testSets;
            X_train_L = mean(RL(train_indices, 2001:end), 1); % 1 x T
            X_test_L = mean(RL(test_indices, 2001:end), 1); % 1 x T

            [N1, T1] = size(X_train_R);
            [N2, T1] = size(X_train_L);

            X_train = [X_train_R'; X_train_L']; % (N1+N2) x T
            Y_train = [repmat("Right", T1, 1) ; ...
                       repmat("Left" , T2, 1)]; % 2T x 1

            %%% metrics w/o/ LDA
            if mean(X_train_R) > mean(X_train_L)
                posLabel = "Right";
            else
                posLabel = "Left";
            end

            scores = X_train;
            [fpr,tpr,~,auROC] = perfcurve(Y_train, X_train, posLabel);
            auroc_incor = [auroc_incor; auROC];

            fpr_list_incor{repeat_idx} = fpr;
            tpr_list_incor{repeat_idx} = tpr;
        end

        selectivity_index(cond_idx).unit(unit_idx).tpr_list_cor = tpr_list_cor;
        selectivity_index(cond_idx).unit(unit_idx).fpr_list_cor = fpr_list_cor;
        selectivity_index(cond_idx).unit(unit_idx).tpr_list_incor = tpr_list_incor;
        selectivity_index(cond_idx).unit(unit_idx).fpr_list_incor = fpr_list_incor;

        selectivity_index(cond_idx).unit(unit_idx).auroc_cor = auroc_cor;
        selectivity_index(cond_idx).unit(unit_idx).auroc_incor = auroc_incor;
    end
    toc

end

save(strcat(param.save_path, '\', 'selective_index_', param.version, '-0toEnd_AllEpochs-avgFR_wo_LDA_repeat', num2str(param.num_repeat), '.mat'), 'selectivity_index', '-v7.3')
load(strcat(param.save_path, '\', 'selective_index_', param.version, '-0toEnd_AllEpochs-avgFR_wo_LDA_repeat', num2str(param.num_repeat), '.mat'))

%% [BASELINE] LDA analysis w/o LDA — R/L label SHUFFLED (chance level) -- trial-level shuffling
clc
close all

time_bin = -param.ITI_time:param.bin_time:param.t_max_truncated-param.bin_time;

param.num_repeat   = 100;
param.ratio_train  = 0.8;
param.ratio_test   = 1 - param.ratio_train;

selectivity_index_shuffle = struct();

for cond_idx = 1:numel(cond_name_list)

    datasets = struct();
    trials   = struct();
    cv_opts  = struct();
    for sess_idx = 1:num_sess_list(cond_idx)
        for trial_idx = 1:numel(trial_types)
            selected_trial_type = trial_types(trial_idx);

            temp_spks = results(cond_idx).sess(sess_idx).(selected_trial_type).allTrials.vals(:, 1:numel(time_bin), :);
            num_trials = size(temp_spks, 3);
            for repeat_idx = 1:param.num_repeat
                rng(repeat_idx)

                perm = randperm(num_trials);
                ntr  = floor(num_trials * param.ratio_train);
                train_indices = perm(1:ntr);
                test_indices  = perm(ntr+1:end);

                trials(sess_idx).(selected_trial_type).repeat(repeat_idx).testSets  = test_indices;
                trials(sess_idx).(selected_trial_type).repeat(repeat_idx).trainSets = train_indices;
            end
        end
    end

    % Find the number of units
    num_units = 0;
    num_units_list = [];
    for sess_idx = 1:num_sess_list(cond_idx)
        num_units = num_units + size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1);
        num_units_list = [num_units_list; size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1)];
    end
    num_units_cumsum2 = cumsum(num_units_list);
    num_units_cumsum1 = [1; num_units_cumsum2(1:end-1)+1];
    num_units_cumsum  = [num_units_cumsum1, num_units_cumsum2];

    tic
    for unit_idx = 1:num_units

        sess_idx = find(unit_idx >= num_units_cumsum(:,1) & unit_idx <= num_units_cumsum(:,2), 1, 'first');
        local_unit_idx = unit_idx - num_units_cumsum(sess_idx,1) + 1;

        % RR / LL / RL / LR trial matrix l(T_full x nTrials)
        RR = squeeze(results(cond_idx).sess(sess_idx).RR.allTrials.vals(local_unit_idx, :, :));
        LL = squeeze(results(cond_idx).sess(sess_idx).LL.allTrials.vals(local_unit_idx, :, :));
        RL = squeeze(results(cond_idx).sess(sess_idx).RL.allTrials.vals(local_unit_idx, :, :));
        LR = squeeze(results(cond_idx).sess(sess_idx).LR.allTrials.vals(local_unit_idx, :, :));

        RR = Perform_smoothing_fr(RR', param);   % (nTrials x T_full)
        LL = Perform_smoothing_fr(LL', param);
        RL = Perform_smoothing_fr(RL', param);
        LR = Perform_smoothing_fr(LR', param);

        RR = (RR - mean(RR(:, 1:ITI_dur_bin), 2));
        LL = (LL - mean(LL(:, 1:ITI_dur_bin), 2));
        RL = (RL - mean(RL(:, 1:ITI_dur_bin), 2));
        LR = (LR - mean(LR(:, 1:ITI_dur_bin), 2));

        auroc_cor   = [];
        auroc_incor = [];
        fpr_list_cor   = {};
        tpr_list_cor   = {};
        fpr_list_incor = {};
        tpr_list_incor = {};

        for repeat_idx = 1:param.num_repeat
            fprintf("[SHUFFLE] Cond: [%d/%d] | Unit: [%d/%d], Repeat: [%d/%d]\n", ...
                cond_idx, numel(cond_name_list), unit_idx, num_units, repeat_idx, param.num_repeat);

            %%% ---- Correct trials (RR vs LL) : label shuffled ----
            tr_R = trials(sess_idx).RR.repeat(repeat_idx).trainSets;
            tr_L = trials(sess_idx).LL.repeat(repeat_idx).trainSets;
            nR = numel(tr_R);  nL = numel(tr_L);

            % trial-level permutation
            pooled = [RR(tr_R, 2001:end); LL(tr_L, 2001:end)];   % (nR+nL) x T
            rng(repeat_idx + (unit_idx-1)*param.num_repeat + 100000);   % seed for shuffling
            shuf = randperm(nR + nL);
            X_train_R = mean(pooled(shuf(1:nR),    :), 1);   % 1 x T  (pseudo-R)
            X_train_L = mean(pooled(shuf(nR+1:end),:), 1);   % 1 x T  (pseudo-L)

            [~, T1] = size(X_train_R);
            [~, T2] = size(X_train_L);
            X_train = [X_train_R'; X_train_L'];                       % (T1+T2) x 1
            Y_train = [repmat("Right", T1, 1); repmat("Left", T2, 1)];

            if mean(X_train_R) > mean(X_train_L)
                posLabel = "Right";
            else
                posLabel = "Left";
            end
            [fpr,tpr,~,auROC] = perfcurve(Y_train, X_train, posLabel);
            auroc_cor = [auroc_cor; auROC];
            fpr_list_cor{repeat_idx} = fpr;
            tpr_list_cor{repeat_idx} = tpr;

            %%% ---- Incorrect trials (LR vs RL) : label shuffled ----
            tr_LR = trials(sess_idx).LR.repeat(repeat_idx).trainSets;
            tr_RL = trials(sess_idx).RL.repeat(repeat_idx).trainSets;
            nLR = numel(tr_LR);  nRL = numel(tr_RL);

            pooled = [LR(tr_LR, 2001:end); RL(tr_RL, 2001:end)];      % (nLR+nRL) x T
            rng(repeat_idx + (unit_idx-1)*param.num_repeat + 200000); % seed for shuffling
            shuf = randperm(nLR + nRL);
            X_train_R = mean(pooled(shuf(1:nLR),    :), 1);
            X_train_L = mean(pooled(shuf(nLR+1:end),:), 1);

            [~, T1] = size(X_train_R);
            [~, T2] = size(X_train_L);
            X_train = [X_train_R'; X_train_L'];
            Y_train = [repmat("Right", T1, 1); repmat("Left", T2, 1)];

            if mean(X_train_R) > mean(X_train_L)
                posLabel = "Right";
            else
                posLabel = "Left";
            end
            [fpr,tpr,~,auROC] = perfcurve(Y_train, X_train, posLabel);
            auroc_incor = [auroc_incor; auROC];
            fpr_list_incor{repeat_idx} = fpr;
            tpr_list_incor{repeat_idx} = tpr;
        end

        selectivity_index_shuffle(cond_idx).unit(unit_idx).tpr_list_cor   = tpr_list_cor;
        selectivity_index_shuffle(cond_idx).unit(unit_idx).fpr_list_cor   = fpr_list_cor;
        selectivity_index_shuffle(cond_idx).unit(unit_idx).tpr_list_incor = tpr_list_incor;
        selectivity_index_shuffle(cond_idx).unit(unit_idx).fpr_list_incor = fpr_list_incor;
        selectivity_index_shuffle(cond_idx).unit(unit_idx).auroc_cor      = auroc_cor;
        selectivity_index_shuffle(cond_idx).unit(unit_idx).auroc_incor    = auroc_incor;
    end
    toc

end

save(strcat(param.save_path, '\', 'selective_index_', param.version, ...
     '-0toEnd_AllEpochs-avgFR_wo_LDA_SHUFFLE_repeat', num2str(param.num_repeat), '.mat'), ...
     'selectivity_index_shuffle', '-v7.3')

load(strcat(param.save_path, '\', 'selective_index_', param.version, ...
     '-0toEnd_AllEpochs-avgFR_wo_LDA_SHUFFLE_repeat', num2str(param.num_repeat), '.mat'))

%% [BASELINE] LDA analysis w/o LDA — Y_train LABEL SHUFFLED (chance level) -- Y_train label shuffling
clc
close all

time_bin = -param.ITI_time:param.bin_time:param.t_max_truncated-param.bin_time;

param.num_repeat   = 100;
param.ratio_train  = 0.8;
param.ratio_test   = 1 - param.ratio_train;

selectivity_index_shuffle = struct();

for cond_idx = 1:numel(cond_name_list)

    datasets = struct();
    trials   = struct();
    cv_opts  = struct();
    for sess_idx = 1:num_sess_list(cond_idx)
        for trial_idx = 1:numel(trial_types)
            selected_trial_type = trial_types(trial_idx);

            temp_spks = results(cond_idx).sess(sess_idx).(selected_trial_type).allTrials.vals(:, 1:numel(time_bin), :);
            num_trials = size(temp_spks, 3);
            for repeat_idx = 1:param.num_repeat
                rng(repeat_idx)   % same split as the original

                perm = randperm(num_trials);
                ntr  = floor(num_trials * param.ratio_train);
                train_indices = perm(1:ntr);
                test_indices  = perm(ntr+1:end);

                trials(sess_idx).(selected_trial_type).repeat(repeat_idx).testSets  = test_indices;
                trials(sess_idx).(selected_trial_type).repeat(repeat_idx).trainSets = train_indices;
            end
        end
    end

    % Find the number of units
    num_units = 0;
    num_units_list = [];
    for sess_idx = 1:num_sess_list(cond_idx)
        num_units = num_units + size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1);
        num_units_list = [num_units_list; size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1)];
    end
    num_units_cumsum2 = cumsum(num_units_list);
    num_units_cumsum1 = [1; num_units_cumsum2(1:end-1)+1];
    num_units_cumsum  = [num_units_cumsum1, num_units_cumsum2];

    tic
    for unit_idx = 1:num_units

        sess_idx = find(unit_idx >= num_units_cumsum(:,1) & unit_idx <= num_units_cumsum(:,2), 1, 'first');
        local_unit_idx = unit_idx - num_units_cumsum(sess_idx,1) + 1;

        % Load RR / LL / RL / LR trial matrices (T_full x nTrials)
        RR = squeeze(results(cond_idx).sess(sess_idx).RR.allTrials.vals(local_unit_idx, :, :));
        LL = squeeze(results(cond_idx).sess(sess_idx).LL.allTrials.vals(local_unit_idx, :, :));
        RL = squeeze(results(cond_idx).sess(sess_idx).RL.allTrials.vals(local_unit_idx, :, :));
        LR = squeeze(results(cond_idx).sess(sess_idx).LR.allTrials.vals(local_unit_idx, :, :));

        RR = Perform_smoothing_fr(RR', param);   % (nTrials x T_full)
        LL = Perform_smoothing_fr(LL', param);
        RL = Perform_smoothing_fr(RL', param);
        LR = Perform_smoothing_fr(LR', param);

        RR = (RR - mean(RR(:, 1:ITI_dur_bin), 2));
        LL = (LL - mean(LL(:, 1:ITI_dur_bin), 2));
        RL = (RL - mean(RL(:, 1:ITI_dur_bin), 2));
        LR = (LR - mean(LR(:, 1:ITI_dur_bin), 2));

        auroc_cor   = [];
        auroc_incor = [];
        fpr_list_cor   = {};
        tpr_list_cor   = {};
        fpr_list_incor = {};
        tpr_list_incor = {};

        for repeat_idx = 1:param.num_repeat
            fprintf("[SHUFFLE-Ylabel] Cond: [%d/%d] | Unit: [%d/%d], Repeat: [%d/%d]\n", ...
                cond_idx, numel(cond_name_list), unit_idx, num_units, repeat_idx, param.num_repeat);

            %%% ---- Correct trials (RR vs LL) : Y_train label shuffled ----
            tr_R = trials(sess_idx).RR.repeat(repeat_idx).trainSets;
            tr_L = trials(sess_idx).LL.repeat(repeat_idx).trainSets;

            X_train_R = mean(RR(tr_R, 2001:end), 1);   % 1 x T  (true R mean)
            X_train_L = mean(LL(tr_L, 2001:end), 1);   % 1 x T  (true L mean)

            [~, T1] = size(X_train_R);
            [~, T2] = size(X_train_L);
            X_train = [X_train_R'; X_train_L'];                       % (T1+T2) x 1
            Y_train = [repmat("Right", T1, 1); repmat("Left", T2, 1)];

            % >>> Randomly permute labels only (group sizes T1, T2 are preserved) <<<
            rng(repeat_idx + (unit_idx-1)*param.num_repeat + 100000);   % seed for shuffling
            Y_train = Y_train(randperm(numel(Y_train)));

            % Reset posLabel based on the shuffled groups (original rule: group with larger mean -> auROC>=0.5)
            if mean(X_train(Y_train=="Right")) > mean(X_train(Y_train=="Left"))
                posLabel = "Right";
            else
                posLabel = "Left";
            end
            [fpr,tpr,~,auROC] = perfcurve(Y_train, X_train, posLabel);
            auroc_cor = [auroc_cor; auROC];
            fpr_list_cor{repeat_idx} = fpr;
            tpr_list_cor{repeat_idx} = tpr;

            %%% ---- Incorrect trials (LR vs RL) : Y_train label shuffled ----
            tr_LR = trials(sess_idx).LR.repeat(repeat_idx).trainSets;
            tr_RL = trials(sess_idx).RL.repeat(repeat_idx).trainSets;

            X_train_R = mean(LR(tr_LR, 2001:end), 1);   % 1 x T
            X_train_L = mean(RL(tr_RL, 2001:end), 1);   % 1 x T

            [~, T1] = size(X_train_R);
            [~, T2] = size(X_train_L);
            X_train = [X_train_R'; X_train_L'];
            Y_train = [repmat("Right", T1, 1); repmat("Left", T2, 1)];

            % >>> Randomly permute labels only <<<
            rng(repeat_idx + (unit_idx-1)*param.num_repeat + 200000);   % seed for shuffling (separate from correct)
            Y_train = Y_train(randperm(numel(Y_train)));

            if mean(X_train(Y_train=="Right")) > mean(X_train(Y_train=="Left"))
                posLabel = "Right";
            else
                posLabel = "Left";
            end
            [fpr,tpr,~,auROC] = perfcurve(Y_train, X_train, posLabel);
            auroc_incor = [auroc_incor; auROC];
            fpr_list_incor{repeat_idx} = fpr;
            tpr_list_incor{repeat_idx} = tpr;
        end

        selectivity_index_shuffle(cond_idx).unit(unit_idx).tpr_list_cor   = tpr_list_cor;
        selectivity_index_shuffle(cond_idx).unit(unit_idx).fpr_list_cor   = fpr_list_cor;
        selectivity_index_shuffle(cond_idx).unit(unit_idx).tpr_list_incor = tpr_list_incor;
        selectivity_index_shuffle(cond_idx).unit(unit_idx).fpr_list_incor = fpr_list_incor;
        selectivity_index_shuffle(cond_idx).unit(unit_idx).auroc_cor      = auroc_cor;
        selectivity_index_shuffle(cond_idx).unit(unit_idx).auroc_incor    = auroc_incor;
    end
    toc

end

save(strcat(param.save_path, '\', 'selective_index_', param.version, ...
     '-0toEnd_AllEpochs-avgFR_wo_LDA_SHUFFLE-Ylabel_repeat', num2str(param.num_repeat), '.mat'), ...
     'selectivity_index_shuffle', '-v7.3')

load(strcat(param.save_path, '\', 'selective_index_', param.version, ...
     '-0toEnd_AllEpochs-avgFR_wo_LDA_SHUFFLE-Ylabel_repeat', num2str(param.num_repeat), '.mat'))

%% Plot for ROC curve (bootstrapping-averaged -> unit-averaged) STD
clc
close all

%%% For all units
for cond_idx = 1 % :numel(cond_name_list)
    fpr_list_cor_all = {};
    fpr_list_incor_all = {};
    tpr_list_cor_all = {};
    tpr_list_incor_all = {};

    fpr_list_cor_units = {};
    fpr_list_incor_units = {};
    tpr_list_cor_units = {};
    tpr_list_incor_units = {};

    num_units = size(selectivity_index(cond_idx).unit, 2);

    for unit_idx = 1:num_units
        fpr_list_cor = selectivity_index(cond_idx).unit(unit_idx).fpr_list_cor;
        tpr_list_cor = selectivity_index(cond_idx).unit(unit_idx).tpr_list_cor;
        [fpr_cor, tpr_cor] = Interpolate_ROC_curve_v0b(fpr_list_cor, tpr_list_cor, param);
        fpr_list_cor_units{unit_idx} = fpr_cor;
        tpr_list_cor_units{unit_idx} = mean(tpr_cor, 1);

        fpr_list_incor = selectivity_index(cond_idx).unit(unit_idx).fpr_list_incor;
        tpr_list_incor = selectivity_index(cond_idx).unit(unit_idx).tpr_list_incor;
        [fpr_incor, tpr_incor] = Interpolate_ROC_curve_v0b(fpr_list_incor, tpr_list_incor, param);
        fpr_list_incor_units{unit_idx} = fpr_incor;
        tpr_list_incor_units{unit_idx} = mean(tpr_incor, 1);
    end

    [fpr_cor, tpr_cor] = Interpolate_ROC_curve_v0b(fpr_list_cor_units, tpr_list_cor_units, param);
    [fpr_incor, tpr_incor] = Interpolate_ROC_curve_v0b(fpr_list_incor_units, tpr_list_incor_units, param);

    % ROC curve
    meanTPR_cor = mean(tpr_cor,   1)';
    semTPR_cor  = std(tpr_cor, [], 1)'; % ./ sqrt(num_units);
    meanTPR_incor = mean(tpr_incor,   1)';
    semTPR_incor  = std(tpr_incor, [], 1)'; % ./ sqrt(num_units);

    fig = figure;
    hold on;
    h = fill([fpr_cor; flipud(fpr_cor)], [meanTPR_cor+semTPR_cor; flipud(meanTPR_cor)], 'g', 'EdgeColor','none');
    set(h, 'FaceAlpha', 0.2);
    h = fill([fpr_incor; flipud(fpr_incor)], [meanTPR_incor; flipud(meanTPR_incor-semTPR_incor)], 'r', 'EdgeColor','none');
    set(h, 'FaceColor', "#ff8c00")
    set(h, 'FaceAlpha', 0.2);

    plot(fpr_cor, meanTPR_cor, 'LineWidth',2, 'color', 'g');
    plot(fpr_incor, meanTPR_incor, 'LineWidth',2, 'color', '#ff8c00');

    % Diagonal (chance)
    plot([0 1],[0 1],'k--');
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    xlim([0, 1])
    ylim([0, 1])
    set(gca, 'FontSize')
    % title(strcat('Average ROC with \pm SEM for clusterID', {' '}, num2str(cluster_idx)));
    hold off;
    print(fig, strcat(fpath, parts_path, '\', cond_name_list(cond_idx), '-ROC-woLDA_woTresol', '-EachUnit_STD'), '-dsvg')
end

%% Plot for ROC curve (bootstrapping-averaged -> unit-averaged) STD  + SHUFFLE baseline
clc
close all

% If the shuffle results are not in memory, load them first (file saved in the first block)
% load(strcat(param.save_path, '\', 'selective_index_', param.version, ...
%      '-0toEnd_AllEpochs-avgFR_wo_LDA_SHUFFLE_repeat', num2str(param.num_repeat), '.mat'))

%%% For all units
for cond_idx = 1 % :numel(cond_name_list)
    fpr_list_cor_all = {};
    fpr_list_incor_all = {};
    tpr_list_cor_all = {};
    tpr_list_incor_all = {};

    fpr_list_cor_units = {};
    fpr_list_incor_units = {};
    tpr_list_cor_units = {};
    tpr_list_incor_units = {};

    num_units = size(selectivity_index(cond_idx).unit, 2);

    for unit_idx = 1:num_units
        fpr_list_cor = selectivity_index(cond_idx).unit(unit_idx).fpr_list_cor;
        tpr_list_cor = selectivity_index(cond_idx).unit(unit_idx).tpr_list_cor;
        [fpr_cor, tpr_cor] = Interpolate_ROC_curve_v0b(fpr_list_cor, tpr_list_cor, param);
        fpr_list_cor_units{unit_idx} = fpr_cor;
        tpr_list_cor_units{unit_idx} = mean(tpr_cor, 1);

        fpr_list_incor = selectivity_index(cond_idx).unit(unit_idx).fpr_list_incor;
        tpr_list_incor = selectivity_index(cond_idx).unit(unit_idx).tpr_list_incor;
        [fpr_incor, tpr_incor] = Interpolate_ROC_curve_v0b(fpr_list_incor, tpr_list_incor, param);
        fpr_list_incor_units{unit_idx} = fpr_incor;
        tpr_list_incor_units{unit_idx} = mean(tpr_incor, 1);
    end

    [fpr_cor, tpr_cor]     = Interpolate_ROC_curve_v0b(fpr_list_cor_units, tpr_list_cor_units, param);
    [fpr_incor, tpr_incor] = Interpolate_ROC_curve_v0b(fpr_list_incor_units, tpr_list_incor_units, param);

    % ROC curve
    meanTPR_cor = mean(tpr_cor,   1)';
    semTPR_cor  = std(tpr_cor, [], 1)'; % ./ sqrt(num_units);
    meanTPR_incor = mean(tpr_incor,   1)';
    semTPR_incor  = std(tpr_incor, [], 1)'; % ./ sqrt(num_units);

    % ===================== SHUFFLE (chance) baseline =====================
    tpr_list_cor_units_shuf   = {};
    fpr_list_cor_units_shuf   = {};
    tpr_list_incor_units_shuf = {};
    fpr_list_incor_units_shuf = {};

    num_units_shuf = size(selectivity_index_shuffle(cond_idx).unit, 2);
    for unit_idx = 1:num_units_shuf
        fpr_list_cor = selectivity_index_shuffle(cond_idx).unit(unit_idx).fpr_list_cor;
        tpr_list_cor = selectivity_index_shuffle(cond_idx).unit(unit_idx).tpr_list_cor;
        [fpr_cor_s, tpr_cor_s] = Interpolate_ROC_curve_v0b(fpr_list_cor, tpr_list_cor, param);
        fpr_list_cor_units_shuf{unit_idx} = fpr_cor_s;
        tpr_list_cor_units_shuf{unit_idx} = mean(tpr_cor_s, 1);

        fpr_list_incor = selectivity_index_shuffle(cond_idx).unit(unit_idx).fpr_list_incor;
        tpr_list_incor = selectivity_index_shuffle(cond_idx).unit(unit_idx).tpr_list_incor;
        [fpr_incor_s, tpr_incor_s] = Interpolate_ROC_curve_v0b(fpr_list_incor, tpr_list_incor, param);
        fpr_list_incor_units_shuf{unit_idx} = fpr_incor_s;
        tpr_list_incor_units_shuf{unit_idx} = mean(tpr_incor_s, 1);
    end

    % Merge cor + incor shuffle curves into one chance band (pooled, since both are at chance)
    fpr_list_shuf_units = [fpr_list_cor_units_shuf, fpr_list_incor_units_shuf];
    tpr_list_shuf_units = [tpr_list_cor_units_shuf, tpr_list_incor_units_shuf];
    [fpr_shuf, tpr_shuf] = Interpolate_ROC_curve_v0b(fpr_list_shuf_units, tpr_list_shuf_units, param);

    meanTPR_shuf = mean(tpr_shuf,   1)';
    semTPR_shuf  = std(tpr_shuf, [], 1)';
    % (For a gray band matched to correct trials only, use just _cor_units_shuf instead of the pooling above)
    % =====================================================================

    fig = figure;
    hold on;

    % --- Shuffle baseline (gray, drawn first so it stays in the background) ---
    h = fill([fpr_shuf; flipud(fpr_shuf)], ...
             [meanTPR_shuf+semTPR_shuf; flipud(meanTPR_shuf-semTPR_shuf)], ...
             [0.5 0.5 0.5], 'EdgeColor','none');
    set(h, 'FaceAlpha', 0.25);
    hg = plot(fpr_shuf, meanTPR_shuf, 'LineWidth', 1.5, 'color', [0.5 0.5 0.5]);

    % --- Correct (green) ---
    h = fill([fpr_cor; flipud(fpr_cor)], [meanTPR_cor+semTPR_cor; flipud(meanTPR_cor)], 'g', 'EdgeColor','none');
    set(h, 'FaceAlpha', 0.2);
    % --- Incorrect (orange) ---
    h = fill([fpr_incor; flipud(fpr_incor)], [meanTPR_incor; flipud(meanTPR_incor-semTPR_incor)], 'r', 'EdgeColor','none');
    set(h, 'FaceColor', "#ff8c00")
    set(h, 'FaceAlpha', 0.2);

    hc = plot(fpr_cor, meanTPR_cor, 'LineWidth',2, 'color', 'g');
    hi = plot(fpr_incor, meanTPR_incor, 'LineWidth',2, 'color', '#ff8c00');

    % Diagonal (theoretical chance)
    plot([0 1],[0 1],'k--');
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    xlim([0, 1])
    ylim([0, 1])
    set(gca, 'FontSize')
    legend([hc, hi, hg], {'Correct','Incorrect','Shuffle'}, 'Location','southeast', 'Box','off');
    % title(strcat('Average ROC with \pm SEM for clusterID', {' '}, num2str(cluster_idx)));
    hold off;
    print(fig, strcat(fpath, parts_path, '\', cond_name_list(cond_idx), '-ROC-woLDA_woTresol', '-EachUnit_STD_withShuffle'), '-dsvg')
end

% %% Plot for selectivity index w/ normalization
% clc
% close all
% 
% edges1 = 0:0.01:1;
% edges2 = 0:0.01:0.7;
% 
% for cond_idx = 1:numel(cond_name_list)
%     auroc_cor_list = [];
% 
%     for unit_idx = 1:size(selectivity_index(cond_idx).unit, 2)
%         auroc_cor_list = [auroc_cor_list; mean(selectivity_index(cond_idx).unit(unit_idx).auroc_cor)];
%     end
% 
%     [counts, edgesOut] = histcounts(auroc_cor_list, edges1);
%     total = sum(counts);
%     pdf_mass = counts / total;
%     cdf_vals = cumsum(pdf_mass);
%     % binCenters = edgesOut(1:end-1) + diff(edgesOut)/2;
%     % stairs(binCenters, cdf_vals, 'LineWidth', 2);
% 
%     fig = figure;
%     h = histogram('BinEdges', edgesOut, 'BinCounts', pdf_mass); % set bar heights directly
%     h.FaceAlpha = 0.3;
%     xline(0.5, 'r:', 'LineWidth', 1.5)
%     xlim([0, 1])
%     ylim([0, 0.12])
%     yticks([0:0.03:0.12])
%     xlabel('auROC')
%     ylabel('Proportion');
%     set(gca, 'FontSize', 15)
%     print(fig, strcat(fpath, parts_path, '\', cond_name_list(cond_idx), '-auROC(norm)_EachUnits_woLDA'), '-dsvg')
% 
%     auroc_cor_list = abs(auroc_cor_list-0.5)*2;
%     [counts, edgesOut] = histcounts(auroc_cor_list, edges1);
%     total = sum(counts);
%     pdf_mass = counts / total;
%     cdf_vals = cumsum(pdf_mass);
%     % binCenters = edgesOut(1:end-1) + diff(edgesOut)/2;
%     % stairs(binCenters, cdf_vals, 'LineWidth', 2);
% 
%     fig = figure;
%     h = histogram('BinEdges', edgesOut, 'BinCounts', pdf_mass); % set bar heights directly
%     h.FaceAlpha = 0.3;
%     xlim([0, 1])
%     ylim([0, 0.12])
%     yticks([0:0.03:0.12])
%     xlabel('Selectivity index (|auROC-0.5|*2)')
%     ylabel('Proportion');
%     set(gca, 'FontSize', 15)
%     print(fig, strcat(fpath, parts_path, '\', cond_name_list(cond_idx), '-auROC(norm)_selectivityIndex_woLDA'), '-dsvg')
% end

%% Plot for selectivity index w/ normalization  + SHUFFLE baseline (gray)
clc
close all

edges1 = 0:0.01:1;
edges2 = 0:0.01:0.7;

for cond_idx = 1 % :numel(cond_name_list)

    % ---------- Data: per-unit mean auROC ----------
    auroc_cor_list = [];
    for unit_idx = 1:size(selectivity_index(cond_idx).unit, 2)
        auroc_cor_list = [auroc_cor_list; mean(selectivity_index(cond_idx).unit(unit_idx).auroc_cor)];
    end

    % ---------- Shuffle: per-unit mean auROC ----------
    auroc_cor_list_shuf = [];
    for unit_idx = 1:size(selectivity_index_shuffle(cond_idx).unit, 2)
        auroc_cor_list_shuf = [auroc_cor_list_shuf; mean(selectivity_index_shuffle(cond_idx).unit(unit_idx).auroc_cor)];
    end

    % ================= Fig 1: auROC distribution =================
    [counts, edgesOut] = histcounts(auroc_cor_list, edges1);
    pdf_mass = counts / sum(counts);

    [counts_s, ~] = histcounts(auroc_cor_list_shuf, edges1);
    pdf_mass_shuf = counts_s / sum(counts_s);

    fig = figure;
    hold on;
    % shuffle (gray, in the background)
    hs = histogram('BinEdges', edgesOut, 'BinCounts', pdf_mass_shuf);
    hs.FaceColor = [0.5 0.5 0.5]; hs.FaceAlpha = 0.3; hs.EdgeColor = 'none';
    % data
    h = histogram('BinEdges', edgesOut, 'BinCounts', pdf_mass);
    h.FaceAlpha = 0.3;
    xline(0.5, 'r:', 'LineWidth', 1.5)
    xlim([0, 1])
    ylim([0, 0.12])
    yticks([0:0.03:0.12])
    xlabel('auROC')
    ylabel('Proportion');
    legend([h hs], {'Data','Shuffle'}, 'Box','off')
    set(gca, 'FontSize', 15)
    hold off;
    print(fig, strcat(fpath, parts_path, '\', cond_name_list(cond_idx), '-auROC(norm)_EachUnits_woLDA_withShuffle'), '-dsvg')

    % ================= Fig 2: Selectivity index =================
    auroc_cor_list      = abs(auroc_cor_list - 0.5) * 2;
    auroc_cor_list_shuf = abs(auroc_cor_list_shuf - 0.5) * 2;

    [counts, edgesOut] = histcounts(auroc_cor_list, edges1);
    pdf_mass = counts / sum(counts);

    [counts_s, ~] = histcounts(auroc_cor_list_shuf, edges1);
    pdf_mass_shuf = counts_s / sum(counts_s);

    fig = figure;
    hold on;
    % shuffle (gray, in the background)
    hs = histogram('BinEdges', edgesOut, 'BinCounts', pdf_mass_shuf);
    hs.FaceColor = [0.5 0.5 0.5]; hs.FaceAlpha = 0.3; hs.EdgeColor = 'none';
    % data
    h = histogram('BinEdges', edgesOut, 'BinCounts', pdf_mass);
    h.FaceAlpha = 0.3;
    xlim([0, 1])
    ylim([0, 0.12])
    yticks([0:0.03:0.12])
    xlabel('Selectivity index (|auROC-0.5|*2)')
    ylabel('Proportion');
    legend([h hs], {'Data','Shuffle'}, 'Box','off')
    set(gca, 'FontSize', 15)
    hold off;
    print(fig, strcat(fpath, parts_path, '\', cond_name_list(cond_idx), '-auROC(norm)_selectivityIndex_woLDA_withShuffle'), '-dsvg')

end

%% LDA analysis w/o/ time resolution - Single-Single (v0e) - all units
time_bin = -param.ITI_time:param.bin_time:param.t_max_truncated-param.bin_time;

% parameters for pseudo population
param.num_repeat = 100;
param.num_data = 100;
param.ratio_train = 0.8;
param.ratio_test = 1 - param.ratio_train;

lda_results = struct();

for cond_idx = 1 % :numel(cond_name_list)

    datasets = struct();
    trials = struct();
    cv_opts = struct();

    for sess_idx = 1:num_sess_list(cond_idx)
        for trial_idx = 1:numel(trial_types)
            selected_trial_type = trial_types(trial_idx);

            % Read sizes only
            vals = results(cond_idx).sess(sess_idx).(selected_trial_type).allTrials.vals;
            num_trials = size(vals, 3);
            num_units  = size(vals, 1);

            % Number of pseudo-population samples
            nTrain = floor(param.num_data * param.ratio_train);
            nTest  = max(1, floor(param.num_data * (1 - param.ratio_train)));  % avoid 0

            for repeat_idx = 1:param.num_repeat
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
                trials(sess_idx).(selected_trial_type).repeat(repeat_idx).trainSets = train_indices_all;
                trials(sess_idx).(selected_trial_type).repeat(repeat_idx).testSets  = test_indices_all;
            end
        end
    end


    % Find the number of units
    num_units = 0;
    num_units_list = [];
    for sess_idx = 1:num_sess_list(cond_idx)
        num_units = num_units + size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1);
        num_units_list = [num_units_list; size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1)];
    end
    num_units_cumsum2 = cumsum(num_units_list);
    num_units_cumsum1 = [1; num_units_cumsum2(1:end-1)+1];
    num_units_cumsum = [num_units_cumsum1, num_units_cumsum2];

    % Make a pseudo population
    % --- Precompute fixed indices/sizes ---
    nTrain = floor(param.num_data * param.ratio_train);
    nTest  = floor(param.num_data * param.ratio_test);
    twin   = 2001:numel(time_bin);          % time window to average over
    Utot   = num_units_cumsum(end, 2);       % total number of units

    for repeat_idx = 1:param.num_repeat
        for trial_idx = 1:numel(trial_types)
            selected_trial_type = trial_types(trial_idx);

            % Output buffers: final features are already time-averaged (1 value/unit), so only 2D is needed
            train_data = zeros(nTrain, Utot, 'single');
            test_data  = zeros(nTest , Utot, 'single');

            % ---------- Do "preprocessing + time averaging" only once per session ----------
            % For each session s: vals [U x T x R] -> smooth -> baseline subtract -> time-avg
            % Result: avg_by_trial{s}: [U_s x R_s] (units x trials)
            avg_by_trial = cell(num_sess_list(cond_idx), 1);
            units_per_sess = zeros(num_sess_list(cond_idx),1);  % cache of unit count per session

            for sess_idx = 1:num_sess_list(cond_idx)
                vals = results(cond_idx).sess(sess_idx).(selected_trial_type).allTrials.vals;   % [U x T x R]
                [U, T, R] = size(vals);
                units_per_sess(sess_idx) = U;

                % Smoothing along the time axis (dim=2)
                vals = movmean(vals, param.smoothing_pts, 2);

                % Subtract the baseline (1:ITI_dur_bin) mean: implicit expansion, no bsxfun
                mu0 = mean(vals(:, 1:ITI_dur_bin, :), 2);     % [U x 1 x R]
                vals = vals - mu0;                             % [U x T x R]

                % Only the mean over 2001:end is used in the end
                vavg = mean(vals(:, twin, :), 2);              % [U x 1 x R]
                avg_by_trial{sess_idx} = reshape(vavg, U, R);  % [U x R]
            end

            % ---------- Build the training set (only selects trials per unit) ----------
            for data_idx = 1:nTrain
                if mod(data_idx, 10) == 0
                    fprintf('Repeat: [%d/%d], Trial: [%d/%d], Data(train): [%d/%d]\n', ...
                        repeat_idx, param.num_repeat, trial_idx, numel(trial_types), data_idx, nTrain);
                end

                row = zeros(1, Utot, 'single');   % row buffer for one sample (= one pseudo-population)
                for sess_idx = 1:num_sess_list(cond_idx)
                    lb = num_units_cumsum(sess_idx, 1);
                    ub = num_units_cumsum(sess_idx, 2);

                    A = avg_by_trial{sess_idx};   % [U_s x R_s]
                    U = units_per_sess(sess_idx);

                    % Trial indices to use for each unit: (1xU) -> (U x 1)
                    idx_per_unit = trials(sess_idx).(selected_trial_type) ...
                        .repeat(repeat_idx).trainSets(data_idx, :).';
                    % Select all at once via linear indexing
                    lin = sub2ind(size(A), (1:U).', idx_per_unit);
                    row(lb:ub) = A(lin);          % place into the [1 x U_s] segment
                end
                train_data(data_idx, :) = row;
            end

            % ---------- Build the test set (same logic) ----------
            for data_idx = 1:nTest
                row = zeros(1, Utot, 'single');
                for sess_idx = 1:num_sess_list(cond_idx)
                    lb = num_units_cumsum(sess_idx, 1);
                    ub = num_units_cumsum(sess_idx, 2);

                    A = avg_by_trial{sess_idx};   % [U_s x R_s]
                    U = units_per_sess(sess_idx);

                    idx_per_unit = trials(sess_idx).(selected_trial_type) ...
                        .repeat(repeat_idx).testSets(data_idx, :).';
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

save(strcat(param.save_path, '\', 'selective_index_', param.version, '-0toEnd_PopUnits-single-single.mat'), 'lda_results', '-v7.3')
load(strcat(param.save_path, '\', 'selective_index_', param.version, '-0toEnd_PopUnits-single-single.mat'))


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
    print(fig, strcat(fpath, parts_path, '\', cond_name_list(cond_idx), '-LDA_wo_Time_resol_AllUnits_single-single'), '-dsvg');
end

%% [BASELINE-POP] LDA w/o time resolution - Single-Single - all units : TRIAL-LEVEL SHUFFLE
time_bin = -param.ITI_time:param.bin_time:param.t_max_truncated-param.bin_time;

param.num_repeat  = 100;
param.num_data    = 100;
param.ratio_train = 0.8;
param.ratio_test  = 1 - param.ratio_train;

lda_results_shuffle = struct();

for cond_idx = 1:numel(cond_name_list)

    datasets = struct();
    trials   = struct();
    cv_opts  = struct();

    % ---- (same as original) generate train/test pool sampling indices ----
    for sess_idx = 1:num_sess_list(cond_idx)
        for trial_idx = 1:numel(trial_types)
            selected_trial_type = trial_types(trial_idx);
            vals = results(cond_idx).sess(sess_idx).(selected_trial_type).allTrials.vals;
            num_trials = size(vals, 3);
            num_units  = size(vals, 1);
            nTrain = floor(param.num_data * param.ratio_train);
            nTest  = max(1, floor(param.num_data * (1 - param.ratio_train)));
            for repeat_idx = 1:param.num_repeat
                rng(repeat_idx);
                if num_trials >= 2
                    perm = randperm(num_trials);
                    nPoolTrain = max(1, min(num_trials-1, round(num_trials * param.ratio_train)));
                    train_pool = perm(1:nPoolTrain);
                    test_pool  = perm(nPoolTrain+1:end);
                else
                    train_pool = 1; test_pool = 1;
                end
                if ~isempty(train_pool)
                    train_indices_all = uint32(train_pool( randi(numel(train_pool), nTrain, num_units) ));
                else
                    train_indices_all = ones(nTrain, num_units, 'uint32');
                end
                if ~isempty(test_pool)
                    test_indices_all = uint32(test_pool( randi(numel(test_pool), nTest, num_units) ));
                else
                    test_indices_all = ones(nTest, num_units, 'uint32');
                end
                trials(sess_idx).(selected_trial_type).repeat(repeat_idx).trainSets = train_indices_all;
                trials(sess_idx).(selected_trial_type).repeat(repeat_idx).testSets  = test_indices_all;
            end
        end
    end

    % ---- number of units ----
    num_units = 0; num_units_list = [];
    for sess_idx = 1:num_sess_list(cond_idx)
        num_units = num_units + size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1);
        num_units_list = [num_units_list; size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1)];
    end
    num_units_cumsum2 = cumsum(num_units_list);
    num_units_cumsum1 = [1; num_units_cumsum2(1:end-1)+1];
    num_units_cumsum  = [num_units_cumsum1, num_units_cumsum2];

    nTrain = floor(param.num_data * param.ratio_train);
    nTest  = floor(param.num_data * param.ratio_test);
    twin   = 2001:numel(time_bin);
    Utot   = num_units_cumsum(end, 2);

    % ---- (same as original) per-session preprocessing + time averaging: precompute all trial types ----
    avg_all        = cell(num_sess_list(cond_idx), 1);
    units_per_sess = zeros(num_sess_list(cond_idx), 1);
    for sess_idx = 1:num_sess_list(cond_idx)
        avg_all{sess_idx} = struct();
        for trial_idx = 1:numel(trial_types)
            tt = trial_types(trial_idx);
            vals = results(cond_idx).sess(sess_idx).(tt).allTrials.vals; % U x T x R
            [U, ~, R] = size(vals);
            units_per_sess(sess_idx) = U;
            vals = movmean(vals, param.smoothing_pts, 2);
            mu0  = mean(vals(:, 1:ITI_dur_bin, :), 2);
            vals = vals - mu0;
            vavg = mean(vals(:, twin, :), 2);
            avg_all{sess_idx}.(tt) = reshape(vavg, U, R); % U x R
        end
    end

    % ---- Per repeat: build pseudo-populations from trial-level shuffled avg matrices ----
    for repeat_idx = 1:param.num_repeat

        % >>> TRIAL-LEVEL SHUFFLE: pool trials from the two conditions and randomly re-split <<<
        avg_shuf = cell(num_sess_list(cond_idx), 1);
        rng(repeat_idx + (cond_idx-1)*param.num_repeat + 500000); % seed for shuffling
        for sess_idx = 1:num_sess_list(cond_idx)
            avg_shuf{sess_idx} = struct();
            % correct pair (RR, LL)
            A1 = avg_all{sess_idx}.RR;  A2 = avg_all{sess_idx}.LL;
            R1 = size(A1,2);            R2 = size(A2,2);
            poolA = [A1, A2];           pA = randperm(R1+R2);
            avg_shuf{sess_idx}.RR = poolA(:, pA(1:R1));
            avg_shuf{sess_idx}.LL = poolA(:, pA(R1+1:end));
            % incorrect pair (LR, RL)
            B1 = avg_all{sess_idx}.LR;  B2 = avg_all{sess_idx}.RL;
            S1 = size(B1,2);            S2 = size(B2,2);
            poolB = [B1, B2];           pB = randperm(S1+S2);
            avg_shuf{sess_idx}.LR = poolB(:, pB(1:S1));
            avg_shuf{sess_idx}.RL = poolB(:, pB(S1+1:end));
        end

        for trial_idx = 1:numel(trial_types)
            tt = trial_types(trial_idx);
            train_data = zeros(nTrain, Utot, 'single');
            test_data  = zeros(nTest , Utot, 'single');

            for data_idx = 1:nTrain
                row = zeros(1, Utot, 'single');
                for sess_idx = 1:num_sess_list(cond_idx)
                    lb = num_units_cumsum(sess_idx,1); ub = num_units_cumsum(sess_idx,2);
                    A  = avg_shuf{sess_idx}.(tt);      U  = units_per_sess(sess_idx);
                    idx_per_unit = trials(sess_idx).(tt).repeat(repeat_idx).trainSets(data_idx, :).';
                    lin = sub2ind(size(A), (1:U).', idx_per_unit);
                    row(lb:ub) = A(lin);
                end
                train_data(data_idx, :) = row;
            end
            for data_idx = 1:nTest
                row = zeros(1, Utot, 'single');
                for sess_idx = 1:num_sess_list(cond_idx)
                    lb = num_units_cumsum(sess_idx,1); ub = num_units_cumsum(sess_idx,2);
                    A  = avg_shuf{sess_idx}.(tt);      U  = units_per_sess(sess_idx);
                    idx_per_unit = trials(sess_idx).(tt).repeat(repeat_idx).testSets(data_idx, :).';
                    lin = sub2ind(size(A), (1:U).', idx_per_unit);
                    row(lb:ub) = A(lin);
                end
                test_data(data_idx, :) = row;
            end

            datasets(repeat_idx).(tt).train_data = train_data;
            datasets(repeat_idx).(tt).test_data  = test_data;
        end
    end

    % ---- (same as original) LDA + ROC ----
    fpr_list_cor = {}; tpr_list_cor = {}; fpr_list_incor = {}; tpr_list_incor = {};
    for repeat_idx = 1:param.num_repeat
        %%% Correct (RR vs LL)
        X_train_R = datasets(repeat_idx).RR.train_data;  X_test_R = datasets(repeat_idx).RR.test_data;
        X_train_L = datasets(repeat_idx).LL.train_data;  X_test_L = datasets(repeat_idx).LL.test_data;
        [R, ~] = size(X_train_R);
        X_train = [X_train_R; X_train_L];
        Y_train = [repmat("Right", R, 1); repmat("Left", R, 1)];
        ldaMdl = fitcdiscr(X_train, Y_train, 'DiscrimType','pseudoLinear', 'Prior','empirical');
        [R, ~] = size(X_test_R);
        X_test = [X_test_R; X_test_L];
        Y_test = [repmat("Right", R, 1); repmat("Left", R, 1)];
        [~, score] = predict(ldaMdl, X_test);
        posLabel = {'Right'}; [~, posIdx] = ismember(posLabel, ldaMdl.ClassNames);
        [fpr,tpr,~,~] = perfcurve(Y_test, score(:,posIdx), posLabel);
        fpr_list_cor{repeat_idx} = fpr;  tpr_list_cor{repeat_idx} = tpr;

        %%% Incorrect (LR vs RL)
        X_train_R = datasets(repeat_idx).LR.train_data;  X_test_R = datasets(repeat_idx).LR.test_data;
        X_train_L = datasets(repeat_idx).RL.train_data;  X_test_L = datasets(repeat_idx).RL.test_data;
        [R, ~] = size(X_train_R);
        X_train = [X_train_R; X_train_L];
        Y_train = [repmat("Right", R, 1); repmat("Left", R, 1)];
        ldaMdl = fitcdiscr(X_train, Y_train, 'DiscrimType','pseudoLinear', 'Prior','empirical');
        [R, ~] = size(X_test_R);
        X_test = [X_test_R; X_test_L];
        Y_test = [repmat("Right", R, 1); repmat("Left", R, 1)];
        [~, score] = predict(ldaMdl, X_test);
        posLabel = {'Right'}; [~, posIdx] = ismember(posLabel, ldaMdl.ClassNames);
        [fpr,tpr,~,~] = perfcurve(Y_test, score(:,posIdx), posLabel);
        fpr_list_incor{repeat_idx} = fpr;  tpr_list_incor{repeat_idx} = tpr;
    end

    lda_results_shuffle(cond_idx).fpr_list_cor   = fpr_list_cor;
    lda_results_shuffle(cond_idx).tpr_list_cor   = tpr_list_cor;
    lda_results_shuffle(cond_idx).fpr_list_incor = fpr_list_incor;
    lda_results_shuffle(cond_idx).tpr_list_incor = tpr_list_incor;
end

save(strcat(param.save_path, '\', 'selective_index_', param.version, '-0toEnd_PopUnits-single-single_SHUFFLE-trial.mat'), 'lda_results_shuffle', '-v7.3')
load(strcat(param.save_path, '\', 'selective_index_', param.version, '-0toEnd_PopUnits-single-single_SHUFFLE-trial.mat'))

%% [BASELINE-POP] LDA w/o time resolution - Single-Single - all units : Y_train LABEL SHUFFLE
time_bin = -param.ITI_time:param.bin_time:param.t_max_truncated-param.bin_time;

param.num_repeat  = 100;
param.num_data    = 100;
param.ratio_train = 0.8;
param.ratio_test  = 1 - param.ratio_train;

lda_results_shuffle = struct();

for cond_idx = 1:numel(cond_name_list)

    datasets = struct();
    trials   = struct();
    cv_opts  = struct();

    % ---- (same as original) generate train/test pool sampling indices ----
    for sess_idx = 1:num_sess_list(cond_idx)
        for trial_idx = 1:numel(trial_types)
            selected_trial_type = trial_types(trial_idx);
            vals = results(cond_idx).sess(sess_idx).(selected_trial_type).allTrials.vals;
            num_trials = size(vals, 3);
            num_units  = size(vals, 1);
            nTrain = floor(param.num_data * param.ratio_train);
            nTest  = max(1, floor(param.num_data * (1 - param.ratio_train)));
            for repeat_idx = 1:param.num_repeat
                rng(repeat_idx);
                if num_trials >= 2
                    perm = randperm(num_trials);
                    nPoolTrain = max(1, min(num_trials-1, round(num_trials * param.ratio_train)));
                    train_pool = perm(1:nPoolTrain);
                    test_pool  = perm(nPoolTrain+1:end);
                else
                    train_pool = 1; test_pool = 1;
                end
                if ~isempty(train_pool)
                    train_indices_all = uint32(train_pool( randi(numel(train_pool), nTrain, num_units) ));
                else
                    train_indices_all = ones(nTrain, num_units, 'uint32');
                end
                if ~isempty(test_pool)
                    test_indices_all = uint32(test_pool( randi(numel(test_pool), nTest, num_units) ));
                else
                    test_indices_all = ones(nTest, num_units, 'uint32');
                end
                trials(sess_idx).(selected_trial_type).repeat(repeat_idx).trainSets = train_indices_all;
                trials(sess_idx).(selected_trial_type).repeat(repeat_idx).testSets  = test_indices_all;
            end
        end
    end

    % ---- number of units ----
    num_units = 0; num_units_list = [];
    for sess_idx = 1:num_sess_list(cond_idx)
        num_units = num_units + size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1);
        num_units_list = [num_units_list; size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 1)];
    end
    num_units_cumsum2 = cumsum(num_units_list);
    num_units_cumsum1 = [1; num_units_cumsum2(1:end-1)+1];
    num_units_cumsum  = [num_units_cumsum1, num_units_cumsum2];

    nTrain = floor(param.num_data * param.ratio_train);
    nTest  = floor(param.num_data * param.ratio_test);
    twin   = 2001:numel(time_bin);
    Utot   = num_units_cumsum(end, 2);

    % ---- (same as original) per-session preprocessing + time averaging ----
    avg_all        = cell(num_sess_list(cond_idx), 1);
    units_per_sess = zeros(num_sess_list(cond_idx), 1);
    for sess_idx = 1:num_sess_list(cond_idx)
        avg_all{sess_idx} = struct();
        for trial_idx = 1:numel(trial_types)
            tt = trial_types(trial_idx);
            vals = results(cond_idx).sess(sess_idx).(tt).allTrials.vals;
            [U, ~, R] = size(vals);
            units_per_sess(sess_idx) = U;
            vals = movmean(vals, param.smoothing_pts, 2);
            mu0  = mean(vals(:, 1:ITI_dur_bin, :), 2);
            vals = vals - mu0;
            vavg = mean(vals(:, twin, :), 2);
            avg_all{sess_idx}.(tt) = reshape(vavg, U, R);
        end
    end

    % ---- (same as original) build pseudo-populations (no shuffling) ----
    for repeat_idx = 1:param.num_repeat
        for trial_idx = 1:numel(trial_types)
            tt = trial_types(trial_idx);
            train_data = zeros(nTrain, Utot, 'single');
            test_data  = zeros(nTest , Utot, 'single');
            for data_idx = 1:nTrain
                row = zeros(1, Utot, 'single');
                for sess_idx = 1:num_sess_list(cond_idx)
                    lb = num_units_cumsum(sess_idx,1); ub = num_units_cumsum(sess_idx,2);
                    A  = avg_all{sess_idx}.(tt);       U  = units_per_sess(sess_idx);
                    idx_per_unit = trials(sess_idx).(tt).repeat(repeat_idx).trainSets(data_idx, :).';
                    lin = sub2ind(size(A), (1:U).', idx_per_unit);
                    row(lb:ub) = A(lin);
                end
                train_data(data_idx, :) = row;
            end
            for data_idx = 1:nTest
                row = zeros(1, Utot, 'single');
                for sess_idx = 1:num_sess_list(cond_idx)
                    lb = num_units_cumsum(sess_idx,1); ub = num_units_cumsum(sess_idx,2);
                    A  = avg_all{sess_idx}.(tt);       U  = units_per_sess(sess_idx);
                    idx_per_unit = trials(sess_idx).(tt).repeat(repeat_idx).testSets(data_idx, :).';
                    lin = sub2ind(size(A), (1:U).', idx_per_unit);
                    row(lb:ub) = A(lin);
                end
                test_data(data_idx, :) = row;
            end
            datasets(repeat_idx).(tt).train_data = train_data;
            datasets(repeat_idx).(tt).test_data  = test_data;
        end
    end

    % ---- LDA + ROC : Y_train LABEL SHUFFLE ----
    fpr_list_cor = {}; tpr_list_cor = {}; fpr_list_incor = {}; tpr_list_incor = {};
    for repeat_idx = 1:param.num_repeat
        %%% Correct (RR vs LL)
        X_train_R = datasets(repeat_idx).RR.train_data;  X_test_R = datasets(repeat_idx).RR.test_data;
        X_train_L = datasets(repeat_idx).LL.train_data;  X_test_L = datasets(repeat_idx).LL.test_data;
        [R, ~] = size(X_train_R);
        X_train = [X_train_R; X_train_L];
        Y_train = [repmat("Right", R, 1); repmat("Left", R, 1)];
        % >>> Randomly permute Y_train labels only (test labels stay true) <<<
        rng(repeat_idx + (cond_idx-1)*param.num_repeat + 700000);
        Y_train = Y_train(randperm(numel(Y_train)));
        ldaMdl = fitcdiscr(X_train, Y_train, 'DiscrimType','pseudoLinear', 'Prior','empirical');
        [R, ~] = size(X_test_R);
        X_test = [X_test_R; X_test_L];
        Y_test = [repmat("Right", R, 1); repmat("Left", R, 1)];
        [~, score] = predict(ldaMdl, X_test);
        posLabel = {'Right'}; [~, posIdx] = ismember(posLabel, ldaMdl.ClassNames);
        [fpr,tpr,~,~] = perfcurve(Y_test, score(:,posIdx), posLabel);
        fpr_list_cor{repeat_idx} = fpr;  tpr_list_cor{repeat_idx} = tpr;

        %%% Incorrect (LR vs RL)
        X_train_R = datasets(repeat_idx).LR.train_data;  X_test_R = datasets(repeat_idx).LR.test_data;
        X_train_L = datasets(repeat_idx).RL.train_data;  X_test_L = datasets(repeat_idx).RL.test_data;
        [R, ~] = size(X_train_R);
        X_train = [X_train_R; X_train_L];
        Y_train = [repmat("Right", R, 1); repmat("Left", R, 1)];
        % >>> Randomly permute Y_train labels only <<<
        rng(repeat_idx + (cond_idx-1)*param.num_repeat + 800000);
        Y_train = Y_train(randperm(numel(Y_train)));
        ldaMdl = fitcdiscr(X_train, Y_train, 'DiscrimType','pseudoLinear', 'Prior','empirical');
        [R, ~] = size(X_test_R);
        X_test = [X_test_R; X_test_L];
        Y_test = [repmat("Right", R, 1); repmat("Left", R, 1)];
        [~, score] = predict(ldaMdl, X_test);
        posLabel = {'Right'}; [~, posIdx] = ismember(posLabel, ldaMdl.ClassNames);
        [fpr,tpr,~,~] = perfcurve(Y_test, score(:,posIdx), posLabel);
        fpr_list_incor{repeat_idx} = fpr;  tpr_list_incor{repeat_idx} = tpr;
    end

    lda_results_shuffle(cond_idx).fpr_list_cor   = fpr_list_cor;
    lda_results_shuffle(cond_idx).tpr_list_cor   = tpr_list_cor;
    lda_results_shuffle(cond_idx).fpr_list_incor = fpr_list_incor;
    lda_results_shuffle(cond_idx).tpr_list_incor = tpr_list_incor;
end

save(strcat(param.save_path, '\', 'selective_index_', param.version, '-0toEnd_PopUnits-single-single_SHUFFLE-Ylabel.mat'), 'lda_results_shuffle', '-v7.3')
load(strcat(param.save_path, '\', 'selective_index_', param.version, '-0toEnd_PopUnits-single-single_SHUFFLE-Ylabel.mat'))

%% Plot for ROC curve (w/o/ temporal resolution) + SHUFFLE baseline (gray)
% Before running: load the desired shuffle .mat so that lda_results_shuffle is in memory
%   (trial version:  ...PopUnits-single-single_SHUFFLE-trial.mat)
%   (Ylabel version: ...PopUnits-single-single_SHUFFLE-Ylabel.mat)
for cond_idx = 1 % :numel(cond_name_list)

    % ---- Data ----
    fpr_list_cor   = lda_results(cond_idx).fpr_list_cor;
    tpr_list_cor   = lda_results(cond_idx).tpr_list_cor;
    fpr_list_incor = lda_results(cond_idx).fpr_list_incor;
    tpr_list_incor = lda_results(cond_idx).tpr_list_incor;

    % ---- Shuffle (cor+incor pooled into one gray chance band) ----
    fpr_list_shuf = [lda_results_shuffle(cond_idx).fpr_list_cor, lda_results_shuffle(cond_idx).fpr_list_incor];
    tpr_list_shuf = [lda_results_shuffle(cond_idx).tpr_list_cor, lda_results_shuffle(cond_idx).tpr_list_incor];
    [fpr_shuf, tpr_shuf] = Interpolate_ROC_curve(fpr_list_shuf, tpr_list_shuf, param);
    meanTPR_shuf = mean(tpr_shuf, 1)';
    semTPR_shuf  = std(tpr_shuf, 0, 1)' ./ sqrt(size(tpr_shuf,1));

    fig = figure;
    hold on;

    % --- Shuffle baseline (gray, drawn first in the background) ---
    h = fill([fpr_shuf; flipud(fpr_shuf)], [meanTPR_shuf+semTPR_shuf; flipud(meanTPR_shuf-semTPR_shuf)], ...
             [0.5 0.5 0.5], 'EdgeColor','none');
    set(h, 'FaceAlpha', 0.25);
    hg = plot(fpr_shuf, meanTPR_shuf, 'LineWidth', 1.5, 'color', [0.5 0.5 0.5]);

    % --- Correct (green) ---
    [fpr_cor, tpr_cor] = Interpolate_ROC_curve(fpr_list_cor, tpr_list_cor, param);
    meanTPR = mean(tpr_cor, 1)'; semTPR = std(tpr_cor, 0, 1)' ./ sqrt(param.num_repeat);
    hc = plot(fpr_cor, meanTPR, 'LineWidth',2, 'color', 'g');
    h = fill([fpr_cor; flipud(fpr_cor)], [meanTPR+semTPR; flipud(meanTPR-semTPR)], 'g', 'EdgeColor','none');
    set(h, 'FaceAlpha', 0.2);

    % --- Incorrect (orange) ---
    [fpr_incor, tpr_incor] = Interpolate_ROC_curve(fpr_list_incor, tpr_list_incor, param);
    meanTPR = mean(tpr_incor, 1)'; semTPR = std(tpr_incor, 0, 1)' ./ sqrt(param.num_repeat);
    hi = plot(fpr_incor, meanTPR, 'LineWidth',2, 'color', '#ff8c00');
    h = fill([fpr_incor; flipud(fpr_incor)], [meanTPR+semTPR; flipud(meanTPR-semTPR)], 'r', 'EdgeColor','none');
    set(h, 'FaceColor', "#ff8c00"); set(h, 'FaceAlpha', 0.2);

    % Diagonal (theoretical chance)
    plot([0 1],[0 1],'k--');
    xlabel('False Positive Rate'); ylabel('True Positive Rate');
    title('Average ROC with \pm SEM');
    legend([hc hi hg], {'Correct','Incorrect','Shuffle'}, 'Location','southeast', 'Box','off');
    hold off; xlim([0 1]); ylim([0 1]); set(gca, 'FontSize', 15)
    print(fig, strcat(fpath, parts_path, '\', cond_name_list(cond_idx), '-LDA_wo_Time_resol_AllUnits_single-single_withShuffle'), '-dsvg');
end