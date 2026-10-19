function [temp_results] = Calculate_all_even_odd_trials_vector(temp_CR, temp_CL, temp_ER, temp_EL, param)

max_pts = param.t_max_truncated*1/param.bin_time; %
ITI_dur = (1:param.ITI_time*param.fs)'; % 2 [sec]
ITI_dur_bin = numel(Calculate_bins_v2(ITI_dur, param.fs, param.bin_time));

temp_CR_all = mean(temp_CR, 3);
temp_CL_all = mean(temp_CL, 3);
temp_ER_all = mean(temp_ER, 3);
temp_EL_all = mean(temp_EL, 3);
temp_CR_all = Perform_smoothing_fr(temp_CR_all, param);
temp_CL_all = Perform_smoothing_fr(temp_CL_all, param);
temp_ER_all = Perform_smoothing_fr(temp_ER_all, param);
temp_EL_all = Perform_smoothing_fr(temp_EL_all, param);
temp_CR_all = (temp_CR_all - mean(temp_CR_all(:, 1:ITI_dur_bin), 2)) ./ std(temp_CR_all(:, 1:ITI_dur_bin), [], 2);
temp_CL_all = (temp_CL_all - mean(temp_CL_all(:, 1:ITI_dur_bin), 2)) ./ std(temp_CL_all(:, 1:ITI_dur_bin), [], 2);
temp_ER_all = (temp_ER_all - mean(temp_ER_all(:, 1:ITI_dur_bin), 2)) ./ std(temp_ER_all(:, 1:ITI_dur_bin), [], 2);
temp_EL_all = (temp_EL_all - mean(temp_EL_all(:, 1:ITI_dur_bin), 2)) ./ std(temp_EL_all(:, 1:ITI_dur_bin), [], 2);

temp_CR_even = mean(temp_CR(:, :, 2:2:end), 3);
temp_CL_even = mean(temp_CL(:, :, 2:2:end), 3);
temp_ER_even = mean(temp_ER(:, :, 2:2:end), 3);
temp_EL_even = mean(temp_EL(:, :, 2:2:end), 3);
temp_CR_even = Perform_smoothing_fr(temp_CR_even, param);
temp_CL_even = Perform_smoothing_fr(temp_CL_even, param);
temp_ER_even = Perform_smoothing_fr(temp_ER_even, param);
temp_EL_even = Perform_smoothing_fr(temp_EL_even, param);
temp_CR_even = (temp_CR_even - mean(temp_CR_even(:, 1:ITI_dur_bin), 2)) ./ std(temp_CR_even(:, 1:ITI_dur_bin), [], 2);
temp_CL_even = (temp_CL_even - mean(temp_CL_even(:, 1:ITI_dur_bin), 2)) ./ std(temp_CL_even(:, 1:ITI_dur_bin), [], 2);
temp_ER_even = (temp_ER_even - mean(temp_ER_even(:, 1:ITI_dur_bin), 2)) ./ std(temp_ER_even(:, 1:ITI_dur_bin), [], 2);
temp_EL_even = (temp_EL_even - mean(temp_EL_even(:, 1:ITI_dur_bin), 2)) ./ std(temp_EL_even(:, 1:ITI_dur_bin), [], 2);

temp_CR_odd = mean(temp_CR(:, :, 1:2:end), 3);
temp_CL_odd = mean(temp_CL(:, :, 1:2:end), 3);
temp_ER_odd = mean(temp_ER(:, :, 1:2:end), 3);
temp_EL_odd = mean(temp_EL(:, :, 1:2:end), 3);
temp_CR_odd = Perform_smoothing_fr(temp_CR_odd, param);
temp_CL_odd = Perform_smoothing_fr(temp_CL_odd, param);
temp_ER_odd = Perform_smoothing_fr(temp_ER_odd, param);
temp_EL_odd = Perform_smoothing_fr(temp_EL_odd, param);
temp_CR_odd = (temp_CR_odd - mean(temp_CR_odd(:, 1:ITI_dur_bin), 2)) ./ std(temp_CR_odd(:, 1:ITI_dur_bin), [], 2);
temp_CL_odd = (temp_CL_odd - mean(temp_CL_odd(:, 1:ITI_dur_bin), 2)) ./ std(temp_CL_odd(:, 1:ITI_dur_bin), [], 2);
temp_ER_odd = (temp_ER_odd - mean(temp_ER_odd(:, 1:ITI_dur_bin), 2)) ./ std(temp_ER_odd(:, 1:ITI_dur_bin), [], 2);
temp_EL_odd = (temp_EL_odd - mean(temp_EL_odd(:, 1:ITI_dur_bin), 2)) ./ std(temp_EL_odd(:, 1:ITI_dur_bin), [], 2);


unstable_unit_indices = find(isnan(sum(temp_CR_all, 2)) | isnan(sum(temp_CR_even, 2)) | isnan(sum(temp_CR_odd, 2)) | ...
                             isnan(sum(temp_CL_all, 2)) | isnan(sum(temp_CL_even, 2)) | isnan(sum(temp_CL_odd, 2)) | ...
                             isnan(sum(temp_ER_all, 2)) | isnan(sum(temp_ER_even, 2)) | isnan(sum(temp_ER_odd, 2)) | ...
                             isnan(sum(temp_EL_all, 2)) | isnan(sum(temp_EL_even, 2)) | isnan(sum(temp_EL_odd, 2)) | ...
                             isinf(sum(temp_CR_all, 2)) | isinf(sum(temp_CR_even, 2)) | isinf(sum(temp_CR_odd, 2)) | ...
                             isinf(sum(temp_CL_all, 2)) | isinf(sum(temp_CL_even, 2)) | isinf(sum(temp_CL_odd, 2)) | ...
                             isinf(sum(temp_ER_all, 2)) | isinf(sum(temp_ER_even, 2)) | isinf(sum(temp_ER_odd, 2)) | ...
                             isinf(sum(temp_EL_all, 2)) | isinf(sum(temp_EL_even, 2)) | isinf(sum(temp_EL_odd, 2)));

stable_unit_indices = 1:size(temp_CR_all, 1);
stable_unit_indices(unstable_unit_indices) = [];

temp_results.CR_all = temp_CR_all(stable_unit_indices, :);
temp_results.CL_all = temp_CL_all(stable_unit_indices, :);
temp_results.ER_all = temp_ER_all(stable_unit_indices, :);
temp_results.EL_all = temp_EL_all(stable_unit_indices, :);

temp_results.CR_even = temp_CR_even(stable_unit_indices, :);
temp_results.CL_even = temp_CL_even(stable_unit_indices, :);
temp_results.ER_even = temp_ER_even(stable_unit_indices, :);
temp_results.EL_even = temp_EL_even(stable_unit_indices, :);

temp_results.CR_odd = temp_CR_odd(stable_unit_indices, :);
temp_results.CL_odd = temp_CL_odd(stable_unit_indices, :);
temp_results.ER_odd = temp_ER_odd(stable_unit_indices, :);
temp_results.EL_odd = temp_EL_odd(stable_unit_indices, :);

end

