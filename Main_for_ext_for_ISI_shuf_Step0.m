clc
clear
close all

%% Hyperparameters to adjust
param = struct();

% add path NPY library for reading Phy
npy_path = "utils\npy-matlab\npy-matlab\";
addpath(npy_path)

addpath("utils\")

save_path = 'C:\Users\sirom\Documents\KimKimEtAl_muscarnic_neuropixels\';
save_path = 'C:\Users\sirom\Documents\KimKimEtAl_muscarnic_neuropixels\';

% File names to read
txt_fname_list = ["\satelliteviewer\YR006 2023-07-08 11-08-14_Before.txt";
                  "\satelliteviewer\YR006 2023-07-11 11-31-05_Atropine.txt";
                  "\satelliteviewer\YR006 2023-07-15 10-55-10_Before_L.txt";
                  "\satelliteviewer\YR006 2023-07-16 11-19-26_Before.txt";
                  "\satelliteviewer\YR006 2023-07-17 10-25-34_Atropine.txt";
                  "\satelliteviewer\YR007 2023-08-20 12-28-35_Before.txt";
                  "\satelliteviewer\YR007 2023-08-22 11-06-59_atropine.txt";
                  "\satelliteviewer\YR007 2023-08-27 11-31-44_Before.txt";
                  "\satelliteviewer\YR007 2023-08-29 11-03-41_Atropine.txt";
                  "\satelliteviewer\YR008 2023-09-17 11-49-37_Before.txt";
                  "\satelliteviewer\YR008 2023-09-19 11-26-17_atropine.txt";
                  "\satelliteviewer\YR008 2023-09-23 10-29-00_Before_L.txt";
                  "\satelliteviewer\YR008 2023-09-25 11-22-48_atropine.txt";
                  "\satelliteviewer\YR006 2023-07-09 11-31-33_saline.txt";
                  "\satelliteviewer\YR007 2023-08-28 10-51-21_saline.txt";
                  "\satelliteviewer\YR008 2023-09-18 11-57-11_saline.txt";];

ks_fname_list = ["\curation_2503\YR006\YR006_230708\imec0_ks2\";
                "\curation_2503\YR006\YR006_230711\imec0_ks2\";
                "\curation_2503\YR006\YR006_230715\imec0_ks2\";
                "\curation_2503\YR006\YR006_230716\imec0_ks2\";
                "\curation_2503\YR006\YR006_230717\imec0_ks2\";
                "\curation_2503\YR007\YR007_230820\imec0_ks2\";
                "\curation_2503\YR007\YR007_230822\imec0_ks2\";
                "\curation_2503\YR007\YR007_230827\imec0_ks2\";
                "\curation_2503\YR007\YR007_230829\imec0_ks2\";
                "\curation_2503\YR008\YR007_230917\imec0_ks2\";
                "\curation_2503\YR008\YR007_230919\imec0_ks2\";
                "\curation_2503\YR008\YR007_230923\imec0_ks2\";
                "\curation_2503\YR008\YR007_230925\imec0_ks2\";
                "\curation_2503\YR006\YR006_230709\imec0_ks2\"; % saline
                "\curation_2503\YR007\YR007_230828\imec0_ks2\";
                "\curation_2503\YR008\YR007_230918\imec0_ks2\"];

spk_ttl_fname_list = ["\ttl\230708_YR006_g0_tcat.nidq.XD_8_0_500.txt";
                    "\ttl\230711_YR006_g0_tcat.nidq.XD_8_0_500.txt";
                    "\ttl\230715_YR006_g0_tcat.nidq.XD_8_0_500.txt";
                    "\ttl\230716_YR006_g0_tcat.nidq.XD_8_0_500.txt";
                    "\ttl\230717_YR006_g0_tcat.nidq.XD_8_0_500.txt";
                    "\ttl\YR007_230820_g0_tcat.nidq.xd_8_0_500.txt";
                    "\ttl\YR007_230822_g0_tcat.nidq.xd_8_0_500.txt";
                    "\ttl\YR007_230827_g0_tcat.nidq.xd_8_0_500.txt";
                    "\ttl\YR007_230829_g0_tcat.nidq.xd_8_0_500.txt";
                    "\ttl\YR007_230917_g0_tcat.nidq.xd_8_0_500.txt";
                    "\ttl\YR007_230919_g0_tcat.nidq.xd_8_0_500.txt";
                    "\ttl\YR007_230923_g1_tcat.nidq.xd_8_0_500.txt";
                    "\ttl\YR007_230925_g0_tcat.nidq.xd_8_0_500.txt";
                    "\ttl\230709_YR006_g0_tcat.nidq.XD_8_0_500.txt"; % saline
                    "\ttl\YR007_230828_g0_tcat.nidq.xd_8_0_500.txt";
                    "\ttl\YR007_230918_g0_tcat.nidq.xd_8_0_500.txt";

                    ];

beh_ttl_fname_list = ["\ttl\230708_YR006_g0_tcat.nidq.XD_8_2_0.txt";
                    "\ttl\230711_YR006_g0_tcat.nidq.XD_8_2_0.txt";
                    "\ttl\230715_YR006_g0_tcat.nidq.XD_8_2_0.txt";
                    "\ttl\230716_YR006_g0_tcat.nidq.XD_8_2_0.txt";
                    "\ttl\230717_YR006_g0_tcat.nidq.XD_8_2_0.txt";
                    "\ttl\YR007_230820_g0_tcat.nidq.xd_8_2_0.txt";
                    "\ttl\YR007_230822_g0_tcat.nidq.xd_8_2_0.txt";
                    "\ttl\YR007_230827_g0_tcat.nidq.xd_8_2_0.txt";
                    "\ttl\YR007_230829_g0_tcat.nidq.xd_8_2_0.txt";
                    "\ttl\YR007_230917_g0_tcat.nidq.xd_8_2_0.txt";
                    "\ttl\YR007_230919_g0_tcat.nidq.xd_8_2_0.txt";
                    "\ttl\YR007_230923_g1_tcat.nidq.xd_8_2_0.txt";
                    "\ttl\YR007_230925_g0_tcat.nidq.xd_8_2_0.txt";
                    "\ttl\230709_YR006_g0_tcat.nidq.XD_8_2_0.txt"; % saline
                    "\ttl\YR007_230828_g0_tcat.nidq.xd_8_2_0.txt";
                    "\ttl\YR007_230918_g0_tcat.nidq.xd_8_2_0.txt";
                    ];

% condition names to analyze
cond_name_list = ["Before"; "Atropine"; "Saline"];

%% Hyperparameters that most do not need to be adjusted
param = struct();
param.fs = 30000; % [Hz]; sampling rate
param.bin_time = 0.001; % [sec]; temporal binning size to calculate firing rateb
param.smoothing_pts = 100; % [idx]; smoothing data points
param.boxcar.a = 1;
param.boxcar.b = ones(1, param.smoothing_pts);
param.ITI_time = 2; % [sec]; Inter-trial interval duration
param.sample_time = 0.1; % [sec]; timing given a sample
param.delay_time = 2; % [sec]
param.gocue_time = param.sample_time + param.delay_time; % [sec]
param.response_dur = 0.5; % [sec]
param.t_max = param.sample_time + param.delay_time + param.response_dur; % [sec]; maximal end time to analyze

param.epochs.stimulus = [0, param.sample_time];
param.epochs.choice = [param.gocue_time-0.6, param.gocue_time];
param.epochs.action = [param.gocue_time, param.gocue_time+0.6];
param.epochs.outcome = [param.gocue_time+2, param.gocue_time+3];

% postulated probe channel locations for each region
param.loc_name = ["DP"; "IL"; "PrL"; "ACC"; "Motor"]; % Region names
param.prob_loc_list = [0, 14; % DP
    15, 66; % IL
    67, 162; % PrL
    163, 228; % ACC
    229, 270]; % Motor

trial_types = ["RR"; "LL"; "RL"; "LR"];

param.all_trials = 1; % save all trial information

param.save_path = save_path;

param.selected_loc_list = [3, 4, 5];

param.num_shuffling = 10;

param.version = "v16fd";

%% Read data and calculate firing vector
for cond_idx = 2 % 1:numel(cond_name_list)
    r_final = struct();

    fprintf('[Condition index: %d/%d] Calculate firing vector\n', cond_idx, numel(cond_name_list));

    cond_name = cond_name_list(cond_idx);
    read_ks_fname_list = ks_fname_list(contains(txt_fname_list, cond_name) | contains(txt_fname_list, lower(cond_name)));
    read_txt_fname_list = txt_fname_list(contains(txt_fname_list, cond_name) | contains(txt_fname_list, lower(cond_name)));
    read_spk_ttl_fname_list = spk_ttl_fname_list(contains(txt_fname_list, cond_name) | contains(txt_fname_list, lower(cond_name)));
    read_beh_ttl_fname_list = beh_ttl_fname_list(contains(txt_fname_list, cond_name) | contains(txt_fname_list, lower(cond_name)));

    num_file = numel(read_ks_fname_list);

    Data = struct();
    SC_trials = struct('RR', struct(), ...
        'LL', struct(), ...
        'RL', struct(), ...
        'LR', struct());

    for file_idx = 1:num_file
        r = struct();

        fprintf('--[File: %d/%d] Read information \n', file_idx, num_file);
        tic

        %%% Read trial information
        txt_fname = strcat(load_path, read_txt_fname_list(file_idx)); % not function; variable name
        ttl_fname = strcat(load_path, read_beh_ttl_fname_list(file_idx)); % not function; variable name
        trials = Read_trial_information_w_TTL_w_millisecScale(txt_fname, ttl_fname, param);

        %%% (Optional) Plot epochs
        Plot_each_trial_epochs(trials, 0); % 1: plot, 0: do not plot

        %%% Read spike sorting results
        ks_fname = strcat(load_path, read_ks_fname_list(file_idx));
        Data = Read_spike_sorting_results_only_good_v2c(ks_fname, Data, file_idx, param);
        Data(file_idx).num_clus = size(Data(file_idx).SU, 2);

        %%% Find correctly/incorrectly lick right and left trials
        SC_trials = Separate_trials_into_cor_incor_w_TTL(trials, SC_trials, file_idx);
        SC_trials = Extract_trials_given_delay(SC_trials, trial_types, param.delay_time, file_idx);

        fprintf('---#(Units): %d\n', size(Data(file_idx).SU, 2))

        toc
    end
end
