function [trials] = Read_trial_information_w_TTL_w_millisecScale(fname, ttl_fname, param)

%%% Input %%%
% fname: .txt file to read
%%% Output %%%
% trials: struct data which contains trial information about time and logs

opt = detectImportOptions(fname);
opt = setvaropts(opt, 'Type', 'string');
opt.DataLines = [1 Inf];
logs = readtable(fname, opt);

opt = detectImportOptions(ttl_fname);
opt = setvaropts(opt, 'Type', 'string');
opt.DataLines = [1 Inf];
ttl = readtable(ttl_fname, opt);


num_logs = size(logs,1);
session_time = table2array(logs(:, 4));
session_name = table2array(logs(:, 3));
session_vals = table2array(logs(:, 5));
% session_time = char(session_time); % YYYY-MM-DD-HH-MM-SS.SSS
% session_time = session_time(:, 9:end); % HH-MM-SS.SSS
% time = zeros(num_logs, 1);
% for i = 1:num_logs
%     time(i, 1) = str2double(session_time(i, 1:2))*60*60 + str2double(session_time(i, 3:4))*60 + ...
%         str2double(session_time(i, 5:6)) + str2double(session_time(i, 7:9)) * 10^(-3);
% end
% session_time = time - time(1);

ref_time_idx = find(session_name=="refX", 1, 'first');

if numel(ref_time_idx) == 0
    ref_time_idx = find(session_name=="refY", 1, 'first');
end

ref_time = session_time(ref_time_idx);
try
    session_time = str2double(session_time) - str2double(ref_time);
catch
    fprintf('here')
end
session_time = (session_time).*10^(-3); % msec to sec
aligned_time = str2double(ttl.Var1);
clear time

%%% Find trial start/end timing
trial_start_idx = [];
for i = 1:num_logs
    if strcmp(session_name(i), 'trialNum')
        trial_start_idx = [trial_start_idx; i];
    end
end
trial_end_idx = [];
for i = 1:num_logs
    if strcmp(session_name(i), 'Away')
        trial_end_idx = [trial_end_idx; i];
    end
end
num_trials = numel(trial_start_idx);

% %%% Align aduino time
% diff_time = session_time(trial_start_idx(1)) - aligned_time(1);
% session_time = session_time - diff_time;

%%% Eliminate erroneous trials
eliminate_idx1 = [];
for idx = 1:num_trials
    temp_logs = session_name(trial_start_idx(idx):trial_end_idx(idx));
    % whether lickOn is located next to lickportsIN
    portsIN_idx = find(strcmp(temp_logs, 'lickportsIN'));
    if strcmp(temp_logs(portsIN_idx+1), 'lickOn') || strcmp(temp_logs(portsIN_idx+1), 'lick2On')
        continue
    else
        eliminate_idx1 = [eliminate_idx1; idx];
    end
end

eliminate_idx2 = [];
pause_idx = find(strcmp(session_name, 'protocolPauseFlag'));
pause_val = session_vals(pause_idx);
if numel(pause_idx) ~= 0
    pause_start_idx = pause_idx(strcmp(pause_val, '1'));
    pause_end_idx = pause_idx(strcmp(pause_val, '0'));
    num_pause = numel(pause_start_idx);
    for idx = 1:num_pause
        eliminate_lb_idx = pause_start_idx(idx);
        eliminate_ub_idx = pause_end_idx(idx);
        eliminate_idx2 = [eliminate_idx2; idx((trial_start_idx >= eliminate_lb_idx) & (trial_end_idx <= eliminate_ub_idx))];
        
        [~, temp_idx2] = min(abs(trial_start_idx - eliminate_lb_idx));
        eliminate_idx2 = [eliminate_idx2; temp_idx2];
    end
end

eliminate_idx = [eliminate_idx1; eliminate_idx2];
eliminate_idx = unique(eliminate_idx);

selected_ttl_idx = 1:num_trials;
selected_ttl_idx(eliminate_idx) = [];
trial_start_idx(eliminate_idx) = [];
trial_end_idx(eliminate_idx) = [];
num_trials = numel(trial_start_idx);

%%% Read start/end duration for each trial
trials = struct('time', [], ...
    'logs', []);
for idx = 1:num_trials
    lb_idx = trial_start_idx(idx);
    ub_idx = trial_end_idx(idx);
    %%% Align aduino time by TTL pulse
    temp_start_time = session_time(lb_idx);
    temp_end_time = session_time(ub_idx);
    temp_ttl_time = aligned_time(selected_ttl_idx(idx));
    diff_time = temp_ttl_time - temp_start_time;
    aligned_session_time = session_time(lb_idx:ub_idx) + diff_time;

    trials(idx).time = [aligned_session_time(2); aligned_session_time(end)]; % 1: trialNum, 2: cue
    trials(idx).logs = [session_time(lb_idx:ub_idx) + diff_time, session_name(lb_idx:ub_idx), session_vals(lb_idx:ub_idx)];
    trials(idx).ttl_time = [aligned_time(selected_ttl_idx(idx)); aligned_time(selected_ttl_idx(idx))];
end

end

