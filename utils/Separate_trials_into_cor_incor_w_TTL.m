function [SC_trials] = Separate_trials_into_cor_incor_w_TTL(trials, SC_trials, file_idx)

%%% Find correctly/incorrectly lick right and left trials
RR_idx = 0;
LL_idx = 0;
RL_idx = 0;
LR_idx = 0;
num_trials = size(trials, 2);
for idx = 1:num_trials
    correct = sum(strcmp(trials(idx).logs(:, 2), 'Correct'));

    if correct == 0 % incorrect trials
        trialType = trials(idx).logs(strcmp(trials(idx).logs(:, 2), 'trialType'), 3); % find log index which has trialType attribute
        if strcmp(trialType, '0') % it is right <RL>
            RL_idx = RL_idx + 1;
            SC_trials(file_idx).RL(RL_idx).time = trials(idx).time;
            SC_trials(file_idx).RL(RL_idx).logs = trials(idx).logs;
            SC_trials(file_idx).RL(RL_idx).ttl_time = trials(idx).ttl_time;
        elseif strcmp(trialType, '1') % it is left <LR>
            LR_idx = LR_idx + 1;
            SC_trials(file_idx).LR(LR_idx).time = trials(idx).time;
            SC_trials(file_idx).LR(LR_idx).logs = trials(idx).logs;
            SC_trials(file_idx).LR(LR_idx).ttl_time = trials(idx).ttl_time;
        end

    else % correct licking
        trialType = trials(idx).logs(strcmp(trials(idx).logs(:, 2), 'trialType'), 3); % find log index which has trialType attribute

        if strcmp(trialType, '0') % it is right
            RR_idx = RR_idx + 1;
            SC_trials(file_idx).RR(RR_idx).time = trials(idx).time;
            SC_trials(file_idx).RR(RR_idx).logs = trials(idx).logs;
            SC_trials(file_idx).RR(RR_idx).ttl_time = trials(idx).ttl_time;

        elseif strcmp(trialType, '1') % it is left
            LL_idx = LL_idx + 1;
            SC_trials(file_idx).LL(LL_idx).time = trials(idx).time;
            SC_trials(file_idx).LL(LL_idx).logs = trials(idx).logs;
            SC_trials(file_idx).LL(LL_idx).ttl_time = trials(idx).ttl_time;

        end
    end
end

end

