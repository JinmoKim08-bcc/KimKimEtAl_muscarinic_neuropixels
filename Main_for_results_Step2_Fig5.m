clc
clear
close all

%% Hyperparameters to adjust
param = struct();

% add path NPY library for reading Phy
npy_path = "utils\npy-matlab\npy-matlab\";
addpath(npy_path)

save_path = 'C:\Users\sirom\Documents\KimKimEtAl_muscarnic_neuropixels\';
parts_path = '\fig5\';

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

%% 

fpath = strcat(param.save_path, '\', param.version, '\');

sep_trials = ["cor"; "CR"; "CL"; "CR-CL"; "incor"; "ER"; "EL"; "ER-EL"; "CR-ER"; "CL-EL"; "Cor+Incor"];

time_bin = -param.ITI_time:param.bin_time:param.t_max_truncated-param.bin_time;
% time_bin = (0:param.bin_time:param.t_max-param.bin_time)';
z_score_bounds = [0, 3];
time_bounds = [-1, time_bin(end)];

stats = struct();

%% Pearson's correlation
clc
% close all
%%% Calculate temporal angle using within-time variability
t_offset = -param.ITI_time; % param.ITI_time-1; % param.ITI_time

if t_offset == -2
    t_lb_idx = 1;
elseif t_offset == -1
    t_lb_idx = 1*1/param.bin_time;
end
t_ub_idx = (param.ITI_time+param.t_max)*1/param.bin_time;
% time_bin = t_offset:param.bin_time:param.t_max_truncated-param.bin_time;
time_bin = -1:param.bin_time:param.t_max_truncated-param.bin_time;

% clc
% % close all

rates = struct();

for cond_idx = 1:numel(cond_name_list)
    CR_all = [];
    CL_all = [];
    ER_all = [];
    EL_all = [];

    CR_even = [];
    CL_even = [];
    ER_even = [];
    EL_even = [];

    CR_odd = [];
    CL_odd = [];
    ER_odd = [];
    EL_odd = [];

    for sess_idx = 1:num_sess_list(cond_idx)
        temp_CR = results(cond_idx).sess(sess_idx).RR.allTrials.vals(:, t_lb_idx:t_ub_idx, :);
        temp_CL = results(cond_idx).sess(sess_idx).LL.allTrials.vals(:, t_lb_idx:t_ub_idx, :);
        temp_ER = results(cond_idx).sess(sess_idx).RL.allTrials.vals(:, t_lb_idx:t_ub_idx, :);
        temp_EL = results(cond_idx).sess(sess_idx).LR.allTrials.vals(:, t_lb_idx:t_ub_idx, :);

        temp_results = Calculate_all_even_odd_trials_vector(temp_CR, temp_CL, temp_ER, temp_EL, param);

        CR_all = [CR_all; temp_results.CR_all];
        CL_all = [CL_all; temp_results.CL_all];
        ER_all = [ER_all; temp_results.ER_all];
        EL_all = [EL_all; temp_results.EL_all];

        CR_even = [CR_even; temp_results.CR_even];
        CL_even = [CL_even; temp_results.CL_even];
        ER_even = [ER_even; temp_results.ER_even];
        EL_even = [EL_even; temp_results.EL_even];

        CR_odd = [CR_odd; temp_results.CR_odd];
        CL_odd = [CL_odd; temp_results.CL_odd];
        ER_odd = [ER_odd; temp_results.ER_odd];
        EL_odd = [EL_odd; temp_results.EL_odd];
    end
    
    CR_all = CR_all(:, 1001:end);
    CL_all = CL_all(:, 1001:end);
    ER_all = ER_all(:, 1001:end);
    EL_all = EL_all(:, 1001:end);

    CR_even = CR_even(:, 1001:end);
    CL_even = CL_even(:, 1001:end);
    ER_even = ER_even(:, 1001:end);
    EL_even = EL_even(:, 1001:end);

    CR_odd = CR_odd(:, 1001:end);
    CL_odd = CL_odd(:, 1001:end);
    ER_odd = ER_odd(:, 1001:end);
    EL_odd = EL_odd(:, 1001:end);
    
    corr_t1_t2_CR = Calculate_population_pearson_corr_fast(CR_all, CR_all);
    corr_t1_t2_CL = Calculate_population_pearson_corr_fast(CL_all, CL_all);
    corr_t1_t2_ER = Calculate_population_pearson_corr_fast(ER_all, ER_all);
    corr_t1_t2_EL = Calculate_population_pearson_corr_fast(EL_all, EL_all);

    r_cor = (corr_t1_t2_CR + corr_t1_t2_CL)/2;
    fig = figure;
    hold on
    imagesc(time_bin, time_bin, r_cor)
    xline(0, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
    xline(param.sample_time, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
    xline(param.gocue_time, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
    yline(0, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
    yline(param.sample_time, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
    yline(param.gocue_time, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
    xlim([-1, time_bin(end)])
    ylim([-1, time_bin(end)])
    xlabel('Time [sec]')
    ylabel('Time [sec]')
    c = colorbar;
    caxis([0, 1])
    c.Label.String = "Pearson's correlation";
    set(gca, 'FontSize', 15)
    print(fig, strcat(param.save_path, parts_path, '\', cond_name_list(cond_idx), '-PearsonsCorr-Correct'), '-dsvg')

    r_incor = (corr_t1_t2_ER + corr_t1_t2_EL)/2;
    fig = figure;
    hold on
    imagesc(time_bin, time_bin, r_incor)
    xline(0, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
    xline(param.sample_time, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
    xline(param.gocue_time, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
    yline(0, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
    yline(param.sample_time, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
    yline(param.gocue_time, ':', 'Color', [1, 1, 1], 'LineWidth', 3)
    xlim([-1, time_bin(end)])
    ylim([-1, time_bin(end)])
    xlabel('Time [sec]')
    ylabel('Time [sec]')
    c = colorbar;
    caxis([0, 1])
    c.Label.String = "Pearson's correlation";
    set(gca, 'FontSize', 15)
    print(fig, strcat(param.save_path, parts_path, '\', cond_name_list(cond_idx), '-PearsonsCorr-Incorrect'), '-dsvg')
end

%% Pearson's correlation - get deltas w/ (A1) averaged diag 
clc
% close all

param.num_repeat = 100;
param.num_units = 100;

%%% Calculate temporal angle using within-time variability
t_offset = -param.ITI_time; % param.ITI_time-1; % param.ITI_time

if t_offset == -2
    t_lb_idx = 1;
elseif t_offset == -1
    t_lb_idx = 1*1/param.bin_time;
end
t_ub_idx = (param.ITI_time+param.t_max)*1/param.bin_time;
% time_bin = t_offset:param.bin_time:param.t_max_truncated-param.bin_time;
time_bin = -1:param.bin_time:param.t_max_truncated-param.bin_time;

clc
close all

rates = struct();
tic
for cond_idx = 1:numel(cond_name_list)
    CR_all = [];
    CL_all = [];
    ER_all = [];
    EL_all = [];

    CR_even = [];
    CL_even = [];
    ER_even = [];
    EL_even = [];

    CR_odd = [];
    CL_odd = [];
    ER_odd = [];
    EL_odd = [];

    for sess_idx = 1:num_sess_list(cond_idx)
        temp_CR = results(cond_idx).sess(sess_idx).RR.allTrials.vals(:, t_lb_idx:t_ub_idx, :);
        temp_CL = results(cond_idx).sess(sess_idx).LL.allTrials.vals(:, t_lb_idx:t_ub_idx, :);
        temp_ER = results(cond_idx).sess(sess_idx).RL.allTrials.vals(:, t_lb_idx:t_ub_idx, :);
        temp_EL = results(cond_idx).sess(sess_idx).LR.allTrials.vals(:, t_lb_idx:t_ub_idx, :);

        temp_results = Calculate_all_even_odd_trials_vector(temp_CR, temp_CL, temp_ER, temp_EL, param);

        CR_all = [CR_all; temp_results.CR_all];
        CL_all = [CL_all; temp_results.CL_all];
        ER_all = [ER_all; temp_results.ER_all];
        EL_all = [EL_all; temp_results.EL_all];

        CR_even = [CR_even; temp_results.CR_even];
        CL_even = [CL_even; temp_results.CL_even];
        ER_even = [ER_even; temp_results.ER_even];
        EL_even = [EL_even; temp_results.EL_even];

        CR_odd = [CR_odd; temp_results.CR_odd];
        CL_odd = [CL_odd; temp_results.CL_odd];
        ER_odd = [ER_odd; temp_results.ER_odd];
        EL_odd = [EL_odd; temp_results.EL_odd];
    end

    CR_all = CR_all(:, 1001:end);
    CL_all = CL_all(:, 1001:end);
    ER_all = ER_all(:, 1001:end);
    EL_all = EL_all(:, 1001:end);

    CR_even = CR_even(:, 1001:end);
    CL_even = CL_even(:, 1001:end);
    ER_even = ER_even(:, 1001:end);
    EL_even = EL_even(:, 1001:end);

    CR_odd = CR_odd(:, 1001:end);
    CL_odd = CL_odd(:, 1001:end);
    ER_odd = ER_odd(:, 1001:end);
    EL_odd = EL_odd(:, 1001:end);

    num_units = size(CR_all, 1);

    for repeat_idx = 1:param.num_repeat
        rng(repeat_idx)
        selected_unit_indices = randperm(num_units, param.num_units);

        fprintf('Pearson corr | cond: [%d/%d], repeat: [%d/%d]\n', cond_idx, numel(cond_name_list), repeat_idx, param.num_repeat)
        corr_t1_t2_CR = Calculate_population_pearson_corr_fast(CR_all(selected_unit_indices, :), CR_all(selected_unit_indices, :));
        corr_t1_t2_CL = Calculate_population_pearson_corr_fast(CL_all(selected_unit_indices, :), CL_all(selected_unit_indices, :));
        corr_t1_t2_ER = Calculate_population_pearson_corr_fast(ER_all(selected_unit_indices, :), ER_all(selected_unit_indices, :));
        corr_t1_t2_EL = Calculate_population_pearson_corr_fast(EL_all(selected_unit_indices, :), EL_all(selected_unit_indices, :));

        r_cor = (corr_t1_t2_CR + corr_t1_t2_CL)/2;
        r_incor = (corr_t1_t2_ER + corr_t1_t2_EL)/2;

        slice_t_lb_idx = find(abs(time_bin - param.sample_time) <= param.bin_time, 1, 'last');
        slice_t_ub_idx = find(abs(time_bin - param.gocue_time) <= param.bin_time, 1, 'first');

        T = slice_t_ub_idx - slice_t_lb_idx + 1;


        %%% Correct trials
        deltas = zeros(T, 1);
        temp_x = r_cor(slice_t_lb_idx:slice_t_ub_idx, slice_t_lb_idx:slice_t_ub_idx);

        for k = 0:-1:-(T-1)
            diag_elements = diag(temp_x, k); % rotate with flipud to get the desired orientation (same orientation as the heatmap plot)
            diag_elements = mean(diag_elements);
            deltas(abs(k)+1) = diag_elements;
        end
        results(cond_idx).repeat(repeat_idx).deltas = deltas;

        %%% Incorrect trials
        deltas = zeros(T, 1);
        temp_x = r_incor(slice_t_lb_idx:slice_t_ub_idx, slice_t_lb_idx:slice_t_ub_idx);

        for k = 0:-1:-(T-1)
            diag_elements = diag(temp_x, k); % rotate with flipud to get the desired orientation (same orientation as the heatmap plot)
            diag_elements = mean(diag_elements);
            deltas(abs(k)+1) = diag_elements;
        end
        results(cond_idx).repeat(repeat_idx).deltas_incor = deltas;
    end    
end
toc

save(strcat(param.save_path, 'pearson_corr_fig5.mat'), 'results', '-v7.3')
load(strcat(param.save_path, 'pearson_corr_fig5.mat'))

%% Pearson's correlation - Delta plot
min_indices_all = [];
for cond_idx = 1:numel(cond_name_list)
    min_indices = zeros(param.num_repeat, 1);
    for repeat_idx = 1:param.num_repeat
        [~, temp_min_idx] = min(results(cond_idx).repeat(repeat_idx).deltas);
        
        min_indices(repeat_idx, 1) = temp_min_idx;

    end

    min_indices_all = [min_indices_all; min_indices];
end

delta_plot_max = median(min_indices_all);

delta_x = 1:delta_plot_max;
delta_x = delta_x';

%%% SEM
fig = figure;
hold on

for cond_idx = 1:numel(cond_name_list)
    delta_all = zeros(delta_plot_max, param.num_repeat);

    for repeat_idx = 1:param.num_repeat
        delta_all(:, repeat_idx) = results(cond_idx).repeat(repeat_idx).deltas(1:delta_plot_max);
    end
    
    delta_mean = mean(delta_all, 2);
    delta_sem = std(delta_all, [], 2) / sqrt(param.num_repeat);
    
    if cond_idx == 1
        h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_sem; flipud(delta_mean-delta_sem)], 'g', 'EdgeColor','none');
        set(h, 'FaceAlpha', 0.2);
        plot(delta_x, delta_mean, 'g', 'LineWidth', 1.5);

    elseif cond_idx == 2
        h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_sem; flipud(delta_mean-delta_sem)], 'g', 'EdgeColor','none');
        set(h, 'FaceColor', "#ff8c00")
        set(h, 'FaceAlpha', 0.2);
        plot(delta_x, delta_mean, 'color', '#ff8c00', 'LineWidth', 1.5);
        
    elseif cond_idx == 3
        h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_sem; flipud(delta_mean-delta_sem)], 'g', 'EdgeColor','none');
        set(h, 'FaceColor', "#7e2f8e")
        set(h, 'FaceAlpha', 0.2);
        plot(delta_x, delta_mean, 'color', '#7e2f8e', 'LineWidth', 1.5);
        
    end
end

% incorrect
delta_all = zeros(delta_plot_max, param.num_repeat);
for repeat_idx = 1:param.num_repeat
    delta_all(:, repeat_idx) = results(1).repeat(repeat_idx).deltas_incor(1:delta_plot_max);
end
delta_mean = mean(delta_all, 2);
delta_sem = std(delta_all, [], 2) / sqrt(param.num_repeat);
h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_sem; flipud(delta_mean-delta_sem)], 'b', 'EdgeColor','none');
set(h, 'FaceColor', [0 0 1.0])
set(h, 'FaceAlpha', 0.2);
plot(delta_x, delta_mean, 'color', 'b', 'LineWidth', 1.5);

% legends
xlabel('\Delta [msec]')
ylabel('Averaged correlation coefficient')
set(gca, 'FontSize', 15)
xlim([0, max(delta_x)])
xticks(0:200:1300)
print(fig, strcat(param.save_path, parts_path, '\', 'Avg corr coeff-delta-SEM error'), '-dsvg');

%%% SEM ---- only for 20 bootstrapping
temp_num_repeat = 20;

fig = figure;
hold on

for cond_idx = 1:numel(cond_name_list)
    delta_all = zeros(delta_plot_max, temp_num_repeat);

    for repeat_idx = 1:temp_num_repeat
        delta_all(:, repeat_idx) = results(cond_idx).repeat(repeat_idx).deltas(1:delta_plot_max);
    end
    
    delta_mean = mean(delta_all, 2);
    delta_sem = std(delta_all, [], 2) / sqrt(temp_num_repeat);
    
    if cond_idx == 1
        h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_sem; flipud(delta_mean-delta_sem)], 'g', 'EdgeColor','none');
        set(h, 'FaceAlpha', 0.2);
        plot(delta_x, delta_mean, 'g', 'LineWidth', 1.5);

    elseif cond_idx == 2
        h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_sem; flipud(delta_mean-delta_sem)], 'g', 'EdgeColor','none');
        set(h, 'FaceColor', "#ff8c00")
        set(h, 'FaceAlpha', 0.2);
        plot(delta_x, delta_mean, 'color', '#ff8c00', 'LineWidth', 1.5);
        
    elseif cond_idx == 3
        h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_sem; flipud(delta_mean-delta_sem)], 'g', 'EdgeColor','none');
        set(h, 'FaceColor', "#7e2f8e")
        set(h, 'FaceAlpha', 0.2);
        plot(delta_x, delta_mean, 'color', '#7e2f8e', 'LineWidth', 1.5);
        
    end
end

% incorrect
delta_all = zeros(delta_plot_max, temp_num_repeat);
for repeat_idx = 1:temp_num_repeat
    delta_all(:, repeat_idx) = results(1).repeat(repeat_idx).deltas_incor(1:delta_plot_max);
end
delta_mean = mean(delta_all, 2);
delta_sem = std(delta_all, [], 2) / sqrt(temp_num_repeat);
h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_sem; flipud(delta_mean-delta_sem)], 'b', 'EdgeColor','none');
set(h, 'FaceColor', [0 0 1.0])
set(h, 'FaceAlpha', 0.2);
plot(delta_x, delta_mean, 'color', 'b', 'LineWidth', 1.5);

% legends
xlabel('\Delta [msec]')
ylabel('Averaged correlation coefficient')
set(gca, 'FontSize', 15)
xlim([0, max(delta_x)])
xticks(0:200:1300)
print(fig, strcat(param.save_path, parts_path, '\', 'Avg corr coeff-delta-SEM error_20bootstrappings'), '-dsvg');


%%% SEM 95%-CI
min_indices_all = [];
for cond_idx = 1:numel(cond_name_list)
    min_indices = zeros(param.num_repeat, 1);
    for repeat_idx = 1:param.num_repeat
        [~, temp_min_idx] = min(results(cond_idx).repeat(repeat_idx).deltas);
        
        min_indices(repeat_idx, 1) = temp_min_idx;
    end

    min_indices_all = [min_indices_all; min_indices];
end

delta_plot_max = median(min_indices_all);

delta_x = 1:delta_plot_max;
delta_x = delta_x';

%%% 95% Confidence Interval (CI)
fig = figure;
hold on

% Settings for computing the t-value for the 95% CI
alpha = 0.05; % 95% confidence level
% Degrees of freedom
df = param.num_repeat - 1; 
% Compute the t-score (two-sided test)
t_score = tinv(1 - alpha/2, df); 

for cond_idx = 1:numel(cond_name_list)
    delta_all = zeros(delta_plot_max, param.num_repeat);

    for repeat_idx = 1:param.num_repeat
        delta_all(:, repeat_idx) = results(cond_idx).repeat(repeat_idx).deltas(1:delta_plot_max);
    end
    
    delta_mean = mean(delta_all, 2);
    
    % Compute SEM
    current_sem = std(delta_all, [], 2) / sqrt(param.num_repeat);
    % Compute 95% CI (SEM * t_score)
    delta_ci = current_sem * t_score;
    
    if cond_idx == 1
        % Draw the error band using delta_ci
        h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_ci; flipud(delta_mean-delta_ci)], 'g', 'EdgeColor','none');
        set(h, 'FaceAlpha', 0.2);
        plot(delta_x, delta_mean, 'g', 'LineWidth', 1.5);

    elseif cond_idx == 2
        h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_ci; flipud(delta_mean-delta_ci)], 'g', 'EdgeColor','none');
        set(h, 'FaceColor', "#ff8c00")
        set(h, 'FaceAlpha', 0.2);
        plot(delta_x, delta_mean, 'color', '#ff8c00', 'LineWidth', 1.5);
        
    elseif cond_idx == 3
        h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_ci; flipud(delta_mean-delta_ci)], 'g', 'EdgeColor','none');
        set(h, 'FaceColor', "#7e2f8e")
        set(h, 'FaceAlpha', 0.2);
        plot(delta_x, delta_mean, 'color', '#7e2f8e', 'LineWidth', 1.5);
        
    end
end

% Apply the same changes to the incorrect trials
delta_all = zeros(delta_plot_max, param.num_repeat);
for repeat_idx = 1:param.num_repeat
    delta_all(:, repeat_idx) = results(1).repeat(repeat_idx).deltas_incor(1:delta_plot_max);
end
delta_mean = mean(delta_all, 2);

% Recompute SEM and CI
current_sem = std(delta_all, [], 2) / sqrt(param.num_repeat);
delta_ci = current_sem * t_score;

h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_ci; flipud(delta_mean-delta_ci)], 'b', 'EdgeColor','none');
set(h, 'FaceColor', [0 0 1.0])
set(h, 'FaceAlpha', 0.2);
plot(delta_x, delta_mean, 'color', 'b', 'LineWidth', 1.5);

% legends
xlabel('\Delta [msec]')
ylabel('Averaged correlation coefficient')
set(gca, 'FontSize', 15)
xlim([0, max(delta_x)])
xticks(0:200:1300)

% Save with "95CI" instead of "SEM" in the file name
print(fig, strcat(param.save_path, parts_path, '\',  'Avg corr coeff-delta-95CI'), '-dsvg');

%%% STD
fig = figure;
hold on

for cond_idx = 1:numel(cond_name_list)
    delta_all = zeros(delta_plot_max, param.num_repeat);

    for repeat_idx = 1:param.num_repeat
        delta_all(:, repeat_idx) = results(cond_idx).repeat(repeat_idx).deltas(1:delta_plot_max);
    end

    delta_mean = mean(delta_all, 2);
    delta_sem = std(delta_all, [], 2); % / sqrt(param.num_repeat);

    if cond_idx == 1
        h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_sem; flipud(delta_mean-delta_sem)], 'g', 'EdgeColor','none');
        set(h, 'FaceAlpha', 0.2);
        plot(delta_x, delta_mean, 'g', 'LineWidth', 1.5);

    elseif cond_idx == 2
        h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_sem; flipud(delta_mean-delta_sem)], 'g', 'EdgeColor','none');
        set(h, 'FaceColor', "#ff8c00")
        set(h, 'FaceAlpha', 0.2);
        plot(delta_x, delta_mean, 'color', '#ff8c00', 'LineWidth', 1.5);

    elseif cond_idx == 3
        h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_sem; flipud(delta_mean-delta_sem)], 'g', 'EdgeColor','none');
        set(h, 'FaceColor', "#7e2f8e")
        set(h, 'FaceAlpha', 0.2);
        plot(delta_x, delta_mean, 'color', '#7e2f8e', 'LineWidth', 1.5);

    end
end

% incorrect
delta_all = zeros(delta_plot_max, param.num_repeat);
for repeat_idx = 1:param.num_repeat
    delta_all(:, repeat_idx) = results(1).repeat(repeat_idx).deltas_incor(1:delta_plot_max);
end
delta_mean = mean(delta_all, 2);
delta_sem = std(delta_all, [], 2); % / sqrt(param.num_repeat);
h = fill([delta_x; flipud(delta_x)], [delta_mean+delta_sem; flipud(delta_mean-delta_sem)], 'b', 'EdgeColor','none');
set(h, 'FaceColor', [0 0 1.0])
set(h, 'FaceAlpha', 0.2);
plot(delta_x, delta_mean, 'color', 'b', 'LineWidth', 1.5);

xlim([0, max(delta_x)])
xticks(0:200:1300)
xlabel('\Delta [ms]')
ylabel('Averaged correlation coefficient')
set(gca, 'FontSize', 15)
print(fig, strcat(param.save_path, parts_path, '\',  'Avg corr coeff-delta-STD error'), '-dsvg');

%% Decay constant (w/ incorrect) - (A1) averaged diag & (B1) exp fitting
% [bestIdx, fitResult, gof, residuals] = find_good_decay_interval(t, y)
clc
close all

max_use_delta = 1000;
min_use_y = 0.3;

for cond_idx = 1:numel(cond_name_list)
    tau = zeros(param.num_repeat, 1);
    R_squared = zeros(param.num_repeat, 1);
    for repeat_idx = 1:param.num_repeat
        [~, temp_min_idx] = min(results(cond_idx).repeat(repeat_idx).deltas);

        input_delta = results(cond_idx).repeat(repeat_idx).deltas(1:temp_min_idx);
        % input_delta = input_delta - min(input_delta)+1e-3;
        input_time = 1:temp_min_idx;

        % [bestIdx, fitResult, gof, residuals] = find_good_decay_interval(input_time, y);

        [t_half, R2] = Estimate_half_life_v0d(input_time, input_delta, min_use_y);

        tau(repeat_idx, 1) = t_half;
        R_squared(repeat_idx, 1) = R2;
    end
    results(cond_idx).tau = tau;
    results(cond_idx).R_squared = R_squared;
end

% incorect
cond_idx = 1;
tau = zeros(param.num_repeat, 1);
R_squared = zeros(param.num_repeat, 1);

for repeat_idx = 1:param.num_repeat
    [~, temp_min_idx] = min(results(cond_idx).repeat(repeat_idx).deltas_incor);

    input_delta = results(cond_idx).repeat(repeat_idx).deltas_incor(1:temp_min_idx);
    % input_delta = input_delta - min(input_delta)+1e-3;
    input_time = 1:temp_min_idx;
    [t_half, R2] = Estimate_half_life_v0d(input_time, input_delta, min_use_y);

    tau(repeat_idx, 1) = t_half;
    R_squared(repeat_idx, 1) = R2;

end
results(cond_idx).tau_incor = tau;
results(cond_idx).R_squared_incor = R_squared;

mean(results(1).R_squared)
mean(results(2).R_squared)
mean(results(3).R_squared)
mean(results(1).R_squared_incor)

%%% SEM
fig = figure;
hold on
for cond_idx = 1:numel(cond_name_list)
    if cond_idx == 1
        bar(1, mean(results(cond_idx).tau), 'FaceColor', 'g')
        errorbar(1, mean(results(cond_idx).tau), std((results(cond_idx).tau))/sqrt(param.num_repeat), ...
                 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 30);
    elseif cond_idx == 2
        bar(4, mean(results(cond_idx).tau), 'FaceColor', '#ff8c00')
        errorbar(4, mean(results(cond_idx).tau), std((results(cond_idx).tau))/sqrt(param.num_repeat), ...
                 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 30);
    elseif cond_idx == 3
        bar(3, mean(results(cond_idx).tau), 'FaceColor', '#7e2f8e')
        errorbar(3, mean(results(cond_idx).tau), std((results(cond_idx).tau))/sqrt(param.num_repeat), ...
                 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 30);
    end
end


bar(2, mean(results(1).tau_incor), 'FaceColor', [0 0 1])
errorbar(2, mean(results(1).tau_incor), std((results(1).tau_incor))/sqrt(param.num_repeat), ...
    'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 30);
xticks(1:4)
xticklabels({'Before', 'Incorrect', 'Saline', 'Atropine'})
ylabel('Decay constant')
set(gca, 'FontSize', 15)

print(fig, strcat(param.save_path, parts_path, '\', 'Decay constant-SEM_A1B1'), '-dsvg');

%%% SEM ---- only for 20 bootstrapping
temp_num_repeat = 20;
fig = figure;
hold on
for cond_idx = 1:numel(cond_name_list)
    if cond_idx == 1
        bar(1, mean(results(cond_idx).tau(1:temp_num_repeat)), 'FaceColor', 'g')
        errorbar(1, mean(results(cond_idx).tau(1:temp_num_repeat)), std((results(cond_idx).tau(1:temp_num_repeat)))/sqrt(temp_num_repeat), ...
                 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 30);
    elseif cond_idx == 2
        bar(4, mean(results(cond_idx).tau(1:temp_num_repeat)), 'FaceColor', '#ff8c00')
        errorbar(4, mean(results(cond_idx).tau(1:temp_num_repeat)), std((results(cond_idx).tau(1:temp_num_repeat)))/sqrt(temp_num_repeat), ...
                 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 30);
    elseif cond_idx == 3
        bar(3, mean(results(cond_idx).tau(1:temp_num_repeat)), 'FaceColor', '#7e2f8e')
        errorbar(3, mean(results(cond_idx).tau(1:temp_num_repeat)), std((results(cond_idx).tau(1:temp_num_repeat)))/sqrt(temp_num_repeat), ...
                 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 30);
    end
end


bar(2, mean(results(1).tau_incor(1:temp_num_repeat)), 'FaceColor', [0 0 1])
errorbar(2, mean(results(1).tau_incor(1:temp_num_repeat)), std((results(1).tau_incor(1:temp_num_repeat)))/sqrt(temp_num_repeat), ...
    'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 30);
xticks(1:4)
xticklabels({'Before', 'Incorrect', 'Saline', 'Atropine'})
ylabel('Decay constant')
set(gca, 'FontSize', 15)

print(fig, strcat(param.save_path, parts_path, '\', 'Decay constant-SEM_A1B1_20bootstrappings'), '-dsvg');



%%% 95CI
alpha = 0.05; % 95% confidence level
% Degrees of freedom
df = param.num_repeat - 1; 
% Compute the t-score (two-sided test)
t_score = tinv(1 - alpha/2, df); 

fig = figure;
hold on
for cond_idx = 1:numel(cond_name_list)
    if cond_idx == 1
        bar(1, mean(results(cond_idx).tau), 'FaceColor', 'g')
        errorbar(1, mean(results(cond_idx).tau), std((results(cond_idx).tau))/sqrt(param.num_repeat)*t_score, ...
                 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 30);
    elseif cond_idx == 2
        bar(4, mean(results(cond_idx).tau), 'FaceColor', '#ff8c00')
        errorbar(4, mean(results(cond_idx).tau), std((results(cond_idx).tau))/sqrt(param.num_repeat)*t_score, ...
                 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 30);
    elseif cond_idx == 3
        bar(3, mean(results(cond_idx).tau), 'FaceColor', '#7e2f8e')
        errorbar(3, mean(results(cond_idx).tau), std((results(cond_idx).tau))/sqrt(param.num_repeat)*t_score, ...
                 'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 30);
    end
end


bar(2, mean(results(1).tau_incor), 'FaceColor', [0 0 1])
errorbar(2, mean(results(1).tau_incor), std((results(1).tau_incor))/sqrt(param.num_repeat)*t_score, ...
    'k', 'linestyle', 'none', 'LineWidth', 1.5, 'CapSize', 30);
xticks(1:4)
xticklabels({'Before', 'Incorrect', 'Saline', 'Atropine'})
ylabel('Decay constant')
set(gca, 'FontSize', 15)

print(fig, strcat(param.save_path, parts_path, '\', 'Decay constant-95CI_A1B1'), '-dsvg');

%%% Statistics
% Extract group data
%%% Data extraction and preprocessing
tau1 = results(1).tau;      % condition 1 (Before)
tau2 = results(2).tau;      % condition 2 (Atropine)
tau3 = results(3).tau;      % condition 3 (Saline)
tau4 = results(1).tau_incor; % condition 1 (Before), incorrect

% Convert all to column vectors
tau1 = tau1(:);
tau2 = tau2(:);
tau3 = tau3(:);
tau4 = tau4(:);

%%% Statistical analysis settings (Mann-Whitney U test & Bonferroni correction)
comparisons = {
    'Before correct vs Before incorrect', tau1, tau4;
    'Before correct vs Saline correct', tau1, tau3;
    'Before correct vs Atropine correct', tau1, tau2;
    'Saline correct vs Atropine correct', tau3, tau2;
};

n_comparisons = size(comparisons, 1);
alpha_raw = 0.05;
alpha_bonf = alpha_raw / n_comparisons; % 0.05 / 4 = 0.0125

fprintf('====================================================================================================\n');
fprintf('Statistical Analysis (Two-sided Mann-Whitney U test)\n');
fprintf('Bonferroni Corrected Alpha: %.5f (Original: %.2f; Comparisons: %d)\n', alpha_bonf, alpha_raw, n_comparisons);
fprintf('====================================================================================================\n');
fprintf('%-35s | %-12s | %-14s | %-12s\n', 'Pair', 'p-value', 'p-value (Bonf)', 'Significant?');
fprintf('----------------------------------------------------------------------------------------------------\n');

%%% Run the analysis
for i = 1:n_comparisons
    pair_name = comparisons{i, 1};
    data_A = comparisons{i, 2};
    data_B = comparisons{i, 3};

    % 1. Apply a two-sided Mann-Whitney U test (Wilcoxon rank sum test) to every comparison
    p_val = ranksum(data_A, data_B, 'tail', 'both');

    % 2. Judge significance using the Bonferroni-corrected p-value
    p_bonf = min(p_val * n_comparisons, 1);
    is_significant = p_bonf < alpha_raw;
    if is_significant
        sig_str = 'Yes (*)';
    else
        sig_str = 'No';
    end

    % 3. Print the raw and corrected p-values
    fprintf('%-35s | %.4e   | %.4e     | %-12s\n', pair_name, p_val, p_bonf, sig_str);
end
fprintf('====================================================================================================\n');
fprintf('* Significant when Bonferroni-adjusted p < %.2f (raw p < %.5f).\n', alpha_raw, alpha_bonf);

%%% Statistics ---- only for 20 bootstrapping
% Extract group data
%%% Data extraction and preprocessing
temp_num_repeat = 20;

tau1 = results(1).tau(1:temp_num_repeat);      % condition 1 (Before)
tau2 = results(2).tau(1:temp_num_repeat);      % condition 2 (Atropine)
tau3 = results(3).tau(1:temp_num_repeat);      % condition 3 (Saline)
tau4 = results(1).tau_incor(1:temp_num_repeat); % condition 1 (Before), incorrect

% Convert all to column vectors
tau1 = tau1(:);
tau2 = tau2(:);
tau3 = tau3(:);
tau4 = tau4(:);

%%% Statistical analysis settings (Mann-Whitney U test & Bonferroni correction)
comparisons = {
    'Before correct vs Before incorrect', tau1, tau4;
    'Before correct vs Saline correct', tau1, tau3;
    'Before correct vs Atropine correct', tau1, tau2;
    'Saline correct vs Atropine correct', tau3, tau2;
};

n_comparisons = size(comparisons, 1);
alpha_raw = 0.05;
alpha_bonf = alpha_raw / n_comparisons; % 0.05 / 4 = 0.0125

fprintf('====================================================================================================\n');
fprintf('Statistical Analysis (Two-sided Mann-Whitney U test)\n');
fprintf('Bonferroni Corrected Alpha: %.5f (Original: %.2f; Comparisons: %d)\n', alpha_bonf, alpha_raw, n_comparisons);
fprintf('====================================================================================================\n');
fprintf('%-35s | %-12s | %-14s | %-12s\n', 'Pair', 'p-value', 'p-value (Bonf)', 'Significant?');
fprintf('----------------------------------------------------------------------------------------------------\n');

%%% Run the analysis
for i = 1:n_comparisons
    pair_name = comparisons{i, 1};
    data_A = comparisons{i, 2};
    data_B = comparisons{i, 3};

    % 1. Apply a two-sided Mann-Whitney U test (Wilcoxon rank sum test) to every comparison
    p_val = ranksum(data_A, data_B, 'tail', 'both');

    % 2. Judge significance using the Bonferroni-corrected p-value
    p_bonf = min(p_val * n_comparisons, 1);
    is_significant = p_bonf < alpha_raw;
    if is_significant
        sig_str = 'Yes (*)';
    else
        sig_str = 'No';
    end

    % 3. Print the raw and corrected p-values
    fprintf('%-35s | %.4e   | %.4e     | %-12s\n', pair_name, p_val, p_bonf, sig_str);
end
fprintf('====================================================================================================\n');
fprintf('* Significant when Bonferroni-adjusted p < %.2f (raw p < %.5f).\n', alpha_raw, alpha_bonf);