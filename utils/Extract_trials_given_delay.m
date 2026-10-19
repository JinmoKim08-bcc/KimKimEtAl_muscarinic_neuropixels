function [SC_trials] = Extract_trials_given_delay(SC_trials, trial_types, delay_time, file_idx)

%%% Calculate delay period (duration) and find trials whose delay time is about 2 sec
delay_period = struct('RR', struct(), ...
    'LL', struct(), ...
    'RL', struct(), ...
    'LR', struct());

for trial_type_idx = 1:numel(trial_types)
    selected_trial_type = trial_types(trial_type_idx);
    num_trials = size(SC_trials(file_idx).(selected_trial_type), 2);

    for trial_idx = 1:num_trials
        temp_logs = SC_trials(file_idx).(selected_trial_type)(trial_idx).logs;

        if trial_type_idx == 1 || trial_type_idx == 3
            cue_idx = strcmp(temp_logs(:, 2), 'cue0');
        else
            cue_idx = strcmp(temp_logs(:, 2), 'cue1');
        end

        end_delay_idx = strcmp(temp_logs(:, 2), 'lickportsIN');
        response_idx = (strcmp(temp_logs(:, 2), 'lick2ON') | strcmp(temp_logs(:, 2), 'lick2On') | strcmp(temp_logs(:, 2), 'lickOn') | strcmp(temp_logs(:, 2), 'lickON'));
        cue_time = str2double(temp_logs(cue_idx,1));
        end_delay_time = str2double(temp_logs(end_delay_idx,1));
        response_start_time = min(str2double(temp_logs(response_idx,1)));
        delay_period(file_idx).(selected_trial_type)(trial_idx).delay_period = end_delay_time-cue_time;
        delay_period(file_idx).(selected_trial_type)(trial_idx).end_delay_time = end_delay_time;
        delay_period(file_idx).(selected_trial_type)(trial_idx).response_start_time = response_start_time;
    end
end

backup_SC_trials = SC_trials;
SC_trials = struct('RR', struct(), ...
    'LL', struct(), ...
    'RL', struct(), ...
    'LR', struct());

for trial_type_idx = 1:numel(trial_types)
    selected_trial_type = trial_types(trial_type_idx);
    temp_period_idx = find(([delay_period(file_idx).(selected_trial_type).delay_period])' >= (delay_time+0.05));
    for trial_idx = 1:numel(temp_period_idx)
        SC_trials(file_idx).(selected_trial_type)(trial_idx).time = backup_SC_trials(file_idx).(selected_trial_type)(temp_period_idx(trial_idx)).time;
        SC_trials(file_idx).(selected_trial_type)(trial_idx).end_delay_time = delay_period(file_idx).(selected_trial_type)(temp_period_idx(trial_idx)).end_delay_time;
        SC_trials(file_idx).(selected_trial_type)(trial_idx).response_start_time = delay_period(file_idx).(selected_trial_type)(temp_period_idx(trial_idx)).response_start_time;
        SC_trials(file_idx).(selected_trial_type)(trial_idx).logs = backup_SC_trials(file_idx).(selected_trial_type)(temp_period_idx(trial_idx)).logs;
        SC_trials(file_idx).(selected_trial_type)(trial_idx).ttl_time = backup_SC_trials(file_idx).(selected_trial_type)(temp_period_idx(trial_idx)).ttl_time;
    end
end

end

