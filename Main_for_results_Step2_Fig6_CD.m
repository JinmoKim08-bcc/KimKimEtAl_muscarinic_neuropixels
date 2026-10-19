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

%% Perform coding direction analysis with aligned axis by CV and alinging to shuffled (w/ baseline subtraction)
%%% Find initial results by coding direction analysis
time_bin = -param.ITI_time:param.bin_time:param.t_max_truncated-param.bin_time;
z_score_bounds = [0, 3];
mode_names = ["stim"; "delay"; "response"; "inagaki_stim"; "inagaki_delay"; "ingaki_response"];
param.cd_opt = 'simple_w_ramping_go';
param.num_folds = 5; % 5-fold CV setting (5-fold is used because some cases contain only a small number of trials)
param.num_bootstrap = 100;
cd_results = struct();

for cond_idx = 1:numel(cond_name_list)
    r_direct = struct();
    wt = struct();
    cd = struct();
    for i_mode = 1:3 % 4:6
        projRR_bootstrp = [];
        projLL_bootstrp = [];
        projRL_bootstrp = [];
        projLR_bootstrp = [];
        projCor_bootstrp = [];
        projIncor_bootstrp = [];
        projRR_bootstrp_shuf = [];
        projLL_bootstrp_shuf = [];
        projRL_bootstrp_shuf = [];
        projLR_bootstrp_shuf = [];
        projCor_bootstrp_shuf = [];
        projIncor_bootstrp_shuf = [];

        % [Changed] Create the partitions for each session in advance, before entering the fold loop
        % This ensures that the 10 folds do not overlap and together cover all of the data.
        session_cv = struct();
        for sess_idx = 1:num_sess_list(cond_idx)
            % Get the number of trials for each condition
            n_CR = size(results(cond_idx).sess(sess_idx).RR.allTrials.vals, 3);
            n_CL = size(results(cond_idx).sess(sess_idx).LL.allTrials.vals, 3);
            n_ER = size(results(cond_idx).sess(sess_idx).RL.allTrials.vals, 3);
            n_EL = size(results(cond_idx).sess(sess_idx).LR.allTrials.vals, 3);

            % Create cvpartition objects (exception handling may be needed if there are fewer data than folds)
            % Here we assume there are enough data and create the KFold partitions
            session_cv(sess_idx).c_CR = cvpartition(n_CR, 'KFold', param.num_folds);
            session_cv(sess_idx).c_CL = cvpartition(n_CL, 'KFold', param.num_folds);
            session_cv(sess_idx).c_ER = cvpartition(n_ER, 'KFold', param.num_folds);
            session_cv(sess_idx).c_EL = cvpartition(n_EL, 'KFold', param.num_folds);
        end

        for fold_idx = 1:param.num_folds
            tic
            fprintf("cond: [%d/%d], i_mode: [%d/%d], fold: [%d/%d]\n", cond_idx, numel(cond_name_list), i_mode, 6, fold_idx, param.num_folds)

            % [Changed] Removed the previous rng(bootstrp_idx)
            % (Already applied when creating the partitions above, so the fixed partitions are used inside each fold)

            CR_train = []; CL_train = []; ER_train = []; EL_train = [];
            CR_test = []; CL_test = []; ER_test = []; EL_test = [];

            for sess_idx = 1:num_sess_list(cond_idx)
                temp_CR = results(cond_idx).sess(sess_idx).RR.allTrials.vals(:, 1:numel(time_bin), :);
                temp_CL = results(cond_idx).sess(sess_idx).LL.allTrials.vals(:, 1:numel(time_bin), :);
                temp_ER = results(cond_idx).sess(sess_idx).RL.allTrials.vals(:, 1:numel(time_bin), :);
                temp_EL = results(cond_idx).sess(sess_idx).LR.allTrials.vals(:, 1:numel(time_bin), :);

                % [Changed] Get the indices for the current fold from the precomputed cvpartition objects
                % training() and test() return logical indices.
                train_idx_CR = training(session_cv(sess_idx).c_CR, fold_idx);
                test_idx_CR = test(session_cv(sess_idx).c_CR, fold_idx);

                train_idx_CL = training(session_cv(sess_idx).c_CL, fold_idx);
                test_idx_CL = test(session_cv(sess_idx).c_CL, fold_idx);

                train_idx_ER = training(session_cv(sess_idx).c_ER, fold_idx);
                test_idx_ER = test(session_cv(sess_idx).c_ER, fold_idx);

                train_idx_EL = training(session_cv(sess_idx).c_EL, fold_idx);
                test_idx_EL = test(session_cv(sess_idx).c_EL, fold_idx);

                % Extract data and compute means (using logical indexing)
                CR_train = [CR_train; mean(temp_CR(:, :, train_idx_CR), 3)];
                CL_train = [CL_train; mean(temp_CL(:, :, train_idx_CL), 3)];
                ER_train = [ER_train; mean(temp_ER(:, :, train_idx_ER), 3)];
                EL_train = [EL_train; mean(temp_EL(:, :, train_idx_EL), 3)];

                CR_test = [CR_test; mean(temp_CR(:, :, test_idx_CR), 3)];
                CL_test = [CL_test; mean(temp_CL(:, :, test_idx_CL), 3)];
                ER_test = [ER_test; mean(temp_ER(:, :, test_idx_ER), 3)];
                EL_test = [EL_test; mean(temp_EL(:, :, test_idx_EL), 3)];
            end

            %%% Find axis
            CR_train_smooth = Perform_smoothing_fr(CR_train, param);
            CL_train_smooth = Perform_smoothing_fr(CL_train, param);
            ER_train_smooth = Perform_smoothing_fr(ER_train, param);
            EL_train_smooth = Perform_smoothing_fr(EL_train, param);

            CR_train_smooth = (CR_train_smooth - mean(CR_train_smooth(:, 1:(param.ITI_time*1/param.bin_time)), 2));
            CL_train_smooth = (CL_train_smooth - mean(CL_train_smooth(:, 1:(param.ITI_time*1/param.bin_time)), 2));
            ER_train_smooth = (ER_train_smooth - mean(ER_train_smooth(:, 1:(param.ITI_time*1/param.bin_time)), 2));
            EL_train_smooth = (EL_train_smooth - mean(EL_train_smooth(:, 1:(param.ITI_time*1/param.bin_time)), 2));

            basis = Perform_CD_analysis_v4b(CR_train_smooth, CL_train_smooth, ER_train_smooth, EL_train_smooth, param);

            %%% Apply projection
            CR_test_smooth = Perform_smoothing_fr(CR_test, param);
            CL_test_smooth = Perform_smoothing_fr(CL_test, param);
            ER_test_smooth = Perform_smoothing_fr(ER_test, param);
            EL_test_smooth = Perform_smoothing_fr(EL_test, param);

            CR_test_smooth = (CR_test_smooth - mean(CR_test_smooth(:, 1:(param.ITI_time*1/param.bin_time)), 2));
            CL_test_smooth = (CL_test_smooth - mean(CL_test_smooth(:, 1:(param.ITI_time*1/param.bin_time)), 2));
            ER_test_smooth = (ER_test_smooth - mean(ER_test_smooth(:, 1:(param.ITI_time*1/param.bin_time)), 2));
            EL_test_smooth = (EL_test_smooth - mean(EL_test_smooth(:, 1:(param.ITI_time*1/param.bin_time)), 2));

            cor = (CR_test_smooth+CL_test_smooth)/2;
            incor = (ER_test_smooth+EL_test_smooth)/2;

            projRR = CR_test_smooth'*basis(:, i_mode);
            projLL = CL_test_smooth'*basis(:, i_mode);
            projRL = ER_test_smooth'*basis(:, i_mode);
            projLR = EL_test_smooth'*basis(:, i_mode);
            projCor = cor'*basis(:, i_mode);
            projIncor = incor'*basis(:, i_mode);

            projRR_bootstrp = [projRR_bootstrp; projRR'];
            projLL_bootstrp = [projLL_bootstrp; projLL'];
            projRL_bootstrp = [projRL_bootstrp; projRL'];
            projLR_bootstrp = [projLR_bootstrp; projLR'];
            projCor_bootstrp = [projCor_bootstrp; projCor'];
            projIncor_bootstrp = [projIncor_bootstrp; projIncor'];

            %%% Find shuffled axis
            for boot_idx = 1:param.num_bootstrap
                rng(boot_idx); % fix the seed (for reproducible partition creation)

                fprintf("---bootstrapping: [%d/%d]\n", boot_idx, param.num_bootstrap)

                % (The logic that shuffles the train set's time axis is unrelated to CV, so it is kept as is)
                shuf_idx = randperm(size(CR_train, 2));
                CR_train_shuf = Perform_smoothing_fr(CR_train(:, shuf_idx), param);
                CL_train_shuf = Perform_smoothing_fr(CL_train(:, shuf_idx), param);
                ER_train_shuf = Perform_smoothing_fr(ER_train(:, shuf_idx), param);
                EL_train_shuf = Perform_smoothing_fr(EL_train(:, shuf_idx), param);

                CR_train_shuf = (CR_train_shuf - mean(CR_train_shuf(:, 1:(param.ITI_time*1/param.bin_time)), 2));
                CL_train_shuf = (CL_train_shuf - mean(CL_train_shuf(:, 1:(param.ITI_time*1/param.bin_time)), 2));
                ER_train_shuf = (ER_train_shuf - mean(ER_train_shuf(:, 1:(param.ITI_time*1/param.bin_time)), 2));
                EL_train_shuf = (EL_train_shuf - mean(EL_train_shuf(:, 1:(param.ITI_time*1/param.bin_time)), 2));

                basis_shuf = Perform_CD_analysis_v4b(CR_train_shuf, CL_train_shuf, ER_train_shuf, EL_train_shuf, param);

                %%% Apply projection (Shuffled)
                CR_test_shuf = Perform_smoothing_fr(CR_test(:, shuf_idx), param);
                CL_test_shuf = Perform_smoothing_fr(CL_test(:, shuf_idx), param);
                ER_test_shuf = Perform_smoothing_fr(ER_test(:, shuf_idx), param);
                EL_test_shuf = Perform_smoothing_fr(EL_test(:, shuf_idx), param);

                CR_test_shuf = (CR_test_shuf - mean(CR_test_shuf(:, 1:(param.ITI_time*1/param.bin_time)), 2));
                CL_test_shuf = (CL_test_shuf - mean(CL_test_shuf(:, 1:(param.ITI_time*1/param.bin_time)), 2));
                ER_test_shuf = (ER_test_shuf - mean(ER_test_shuf(:, 1:(param.ITI_time*1/param.bin_time)), 2));
                EL_test_shuf = (EL_test_shuf - mean(EL_test_shuf(:, 1:(param.ITI_time*1/param.bin_time)), 2));

                cor = (CR_test_shuf+CL_test_shuf)/2;
                incor = (ER_test_shuf+EL_test_shuf)/2;

                projRR = CR_test_shuf'*basis_shuf(:, i_mode);
                projLL = CL_test_shuf'*basis_shuf(:, i_mode);
                projRL = ER_test_shuf'*basis_shuf(:, i_mode);
                projLR = EL_test_shuf'*basis_shuf(:, i_mode);
                projCor = cor'*basis_shuf(:, i_mode);
                projIncor = incor'*basis_shuf(:, i_mode);

                projRR_bootstrp_shuf = [projRR_bootstrp_shuf; projRR'];
                projLL_bootstrp_shuf = [projLL_bootstrp_shuf; projLL'];
                projRL_bootstrp_shuf = [projRL_bootstrp_shuf; projRL'];
                projLR_bootstrp_shuf = [projLR_bootstrp_shuf; projLR'];
                projCor_bootstrp_shuf = [projCor_bootstrp_shuf; projCor'];
                projIncor_bootstrp_shuf = [projIncor_bootstrp_shuf; projIncor'];
            end
            toc
        end
        cd_results(cond_idx).modes(i_mode).projRR_bootstrp = projRR_bootstrp;
        cd_results(cond_idx).modes(i_mode).projLL_bootstrp = projLL_bootstrp;
        cd_results(cond_idx).modes(i_mode).projRL_bootstrp = projRL_bootstrp;
        cd_results(cond_idx).modes(i_mode).projLR_bootstrp = projLR_bootstrp;

        cd_results(cond_idx).modes(i_mode).projRR_bootstrp_shuf = projRR_bootstrp_shuf;
        cd_results(cond_idx).modes(i_mode).projLL_bootstrp_shuf = projLL_bootstrp_shuf;
        cd_results(cond_idx).modes(i_mode).projRL_bootstrp_shuf = projRL_bootstrp_shuf;
        cd_results(cond_idx).modes(i_mode).projLR_bootstrp_shuf = projLR_bootstrp_shuf;
    end
end

save(strcat(param.save_path, 'cd_results_v16fd_CV5_vcU.mat'), 'cd_results', '-v7.3')
load(strcat(param.save_path, 'cd_results_v16fd_CV5_vcU.mat'))

%% Step 2. Perform aligning to subspace
subspace_name = ["projRR_bootstrp";
                "projRL_bootstrp";
                "projLR_bootstrp";
                "projLL_bootstrp"];

subspace_name_shuf = ["projRR_bootstrp_shuf";
                      "projRL_bootstrp_shuf";
                      "projLR_bootstrp_shuf";
                      "projLL_bootstrp_shuf"];

for cond_idx = 1:numel(cond_name_list)
    for i_mode = 1:3 %4:6
        for subspace_idx = 1:numel(subspace_name)
            selected_subspace = subspace_name(subspace_idx);
            selected_subspace_shuf = subspace_name_shuf(subspace_idx);

            z_vals = [];

            for fold_idx = 1:param.num_folds
                temp_L_A = cd_results(cond_idx).modes(i_mode).(selected_subspace)(fold_idx, :);

                lb_idx = (fold_idx-1)*param.num_bootstrap + 1;
                ub_idx = fold_idx*param.num_bootstrap;
                temp_L_B = cd_results(cond_idx).modes(i_mode).(selected_subspace_shuf)(lb_idx:ub_idx, :);

                mu_shuf = mean(temp_L_B, 1);
                std_shuf = std(temp_L_B, [], 1);

                z_vals = [z_vals; (temp_L_A - mu_shuf) ./ std_shuf];
            end
            cd_results(cond_idx).modes(i_mode).aligned.(selected_subspace) = z_vals;
        end
    end
end

%% Plot
clc
close all

for cond_idx = 1:3
    for i_mode = 1:3 % 4:6
        projRR_bootstrp = cd_results(cond_idx).modes(i_mode).aligned.projRR_bootstrp;
        projRL_bootstrp = cd_results(cond_idx).modes(i_mode).aligned.projRL_bootstrp;
        projLR_bootstrp = cd_results(cond_idx).modes(i_mode).aligned.projLR_bootstrp;
        projLL_bootstrp = cd_results(cond_idx).modes(i_mode).aligned.projLL_bootstrp;

        fig = figure;
        hold on
        proj_upper = mean(projRR_bootstrp, 1) + std(projRR_bootstrp, [], 1)/sqrt(param.num_bootstrap);
        proj_lower = mean(projRR_bootstrp, 1) - std(projRR_bootstrp, [], 1)/sqrt(param.num_bootstrap);
        h = fill([time_bin'; flipud(time_bin')], [proj_upper'; flipud(proj_lower')], ...
            [0.5, 0.5, 1], 'LineStyle', ':', 'LineWidth', 0.05, 'HandleVisibility', 'off');
        set(h, 'FaceAlpha', 0.3);

        %         proj_upper = mean(projRL_bootstrp, 1) + std(projRL_bootstrp, [], 1)/sqrt(param.num_bootstrap);
        %         proj_lower = mean(projRL_bootstrp, 1) - std(projRL_bootstrp, [], 1)/sqrt(param.num_bootstrap);
        %         h = fill([time_bin'; flipud(time_bin')], [proj_upper'; flipud(proj_lower')], ...
        %             [0.5, 1, 0.5], 'LineStyle', ':', 'LineWidth', 0.05, 'HandleVisibility', 'off');
        %         set(h, 'FaceAlpha', 0.3);

        proj_upper = mean(projLL_bootstrp, 1) + std(projLL_bootstrp, [], 1)/sqrt(param.num_bootstrap);
        proj_lower = mean(projLL_bootstrp, 1) - std(projLL_bootstrp, [], 1)/sqrt(param.num_bootstrap);
        h = fill([time_bin'; flipud(time_bin')], [proj_upper'; flipud(proj_lower')], ...
            [1, 0.5, 0.5], 'LineStyle', ':', 'LineWidth', 0.05, 'HandleVisibility', 'off');
        set(h, 'FaceAlpha', 0.3);

        %         proj_upper = mean(projLR_bootstrp, 1) + std(projLR_bootstrp, [], 1)/sqrt(param.num_bootstrap);
        %         proj_lower = mean(projLR_bootstrp, 1) - std(projLR_bootstrp, [], 1)/sqrt(param.num_bootstrap);
        %         h = fill([time_bin'; flipud(time_bin')], [proj_upper'; flipud(proj_lower')], ...
        %             [1, 0.7745, 0.5], 'LineStyle', ':', 'LineWidth', 0.25, 'HandleVisibility', 'off');
        %         set(h, 'FaceAlpha', 0.3);

        plot(time_bin, mean(projRR_bootstrp, 1), 'b', 'LineWidth', 2)
        %         plot(time_bin, mean(projRL_bootstrp, 1), 'color', [0, 0.8, 0], 'LineWidth', 2)
        plot(time_bin, mean(projLL_bootstrp, 1), 'r', 'LineWidth', 2)
        %         plot(time_bin, mean(projLR_bootstrp, 1), 'color', '#ff8c00', 'LineWidth', 2)
        xline(0, ':', 'Color', [0, 0, 0], 'LineWidth', 3)
        xline(param.sample_time, ':', 'Color', [0, 0, 0], 'LineWidth', 3)
        xline(param.gocue_time, ':', 'Color', [0, 0, 0], 'LineWidth', 3)
        xlim(time_bounds)
        xlabel('Time [sec]')
        ylabel('Projection [a.u.]')
        if i_mode == 1 || i_mode == 4
            % %             legend('RR', 'RL', 'LL', 'LR')
            ylim([-50, 50])
        elseif i_mode == 2 || i_mode == 5
            ylim([-100, 50])
        elseif i_mode == 3 || i_mode == 6
            ylim([-50, 100])
        end
        set(gca, 'FontSize', 15)

        % print(fig, strcat(param.save_path, param.version, '\CD_', param.cd_opt, '_mode', mode_names(i_mode), '_cond', num2str(cond_idx), '_CCA_baseline'), '-dsvg')
    end
end

%% correct trials - distance
clc
close all

for i_mode = 2 % 1:3
    projRR_bootstrp = cd_results(1).modes(i_mode).aligned.projRR_bootstrp;
    projLL_bootstrp = cd_results(1).modes(i_mode).aligned.projLL_bootstrp;
    projDist_before = abs(projRR_bootstrp - projLL_bootstrp);

    projRR_bootstrp = cd_results(2).modes(i_mode).aligned.projRR_bootstrp;
    projLL_bootstrp = cd_results(2).modes(i_mode).aligned.projLL_bootstrp;
    projDist_atropine = abs(projRR_bootstrp - projLL_bootstrp);

    projRR_bootstrp = cd_results(3).modes(i_mode).aligned.projRR_bootstrp;
    projLL_bootstrp = cd_results(3).modes(i_mode).aligned.projLL_bootstrp;
    projDist_saline = abs(projRR_bootstrp - projLL_bootstrp);

    fig = figure;
    hold on
    proj_upper = mean(projDist_before, 1) + std(projDist_before, [], 1)/sqrt(param.num_bootstrap);
    proj_lower = mean(projDist_before, 1) - std(projDist_before, [], 1)/sqrt(param.num_bootstrap);
    h = fill([time_bin'; flipud(time_bin')], [proj_upper'; flipud(proj_lower')], ...
        'g', 'LineStyle', ':', 'LineWidth', 0.05, 'HandleVisibility', 'off');
    set(h, 'FaceAlpha', 0.3);

    proj_upper = mean(projDist_atropine, 1) + std(projDist_atropine, [], 1)/sqrt(param.num_bootstrap);
    proj_lower = mean(projDist_atropine, 1) - std(projDist_atropine, [], 1)/sqrt(param.num_bootstrap);
    h = fill([time_bin'; flipud(time_bin')], [proj_upper'; flipud(proj_lower')], ...
        [250, 144, 0]/255, 'LineStyle', ':', 'LineWidth', 0.05, 'HandleVisibility', 'off');
    set(h, 'FaceAlpha', 0.3);

    proj_upper = mean(projDist_saline, 1) + std(projDist_saline, [], 1)/sqrt(param.num_bootstrap);
    proj_lower = mean(projDist_saline, 1) - std(projDist_saline, [], 1)/sqrt(param.num_bootstrap);
    h = fill([time_bin'; flipud(time_bin')], [proj_upper'; flipud(proj_lower')], ...
        [0.494, 0.184, 0.557], 'LineStyle', ':', 'LineWidth', 0.05, 'HandleVisibility', 'off');
    set(h, 'FaceAlpha', 0.3);

    plot(time_bin, mean(projDist_before, 1), 'g', 'LineWidth', 2)
    plot(time_bin, mean(projDist_atropine, 1), 'color', [255, 144, 0]/255, 'LineWidth', 2)
    plot(time_bin, mean(projDist_saline, 1), 'color', [0.494, 0.184, 0.557], 'LineWidth', 2)

    xline(0, ':', 'Color', [0, 0, 0], 'LineWidth', 3)
    xline(param.sample_time, ':', 'Color', [0, 0, 0], 'LineWidth', 3)
    xline(param.gocue_time, ':', 'Color', [0, 0, 0], 'LineWidth', 3)
    xlim(time_bounds)
    xlabel('Time [sec]')
    ylabel('Projection [a.u.]')
    if i_mode == 1 || i_mode == 4
        % %             legend('RR', 'RL', 'LL', 'LR')
        ylim([-5, 10])
    elseif i_mode == 2 || i_mode == 5
        ylim([-10, 20])
    elseif i_mode == 3 || i_mode == 6
        ylim([-10, 20])
    end
    set(gca, 'FontSize', 15)
    grid on

    print(fig, strcat(param.save_path, parts_path, '\CD_Distance_', param.cd_opt, '_mode', mode_names(i_mode), '_cond', num2str(cond_idx), '_CCA_baseline'), '-dsvg')
end
