function summary_table = Calculate_pooled_summary(value_table, param, value_field, prep_names)

summary_table = table();

for temp_prep_idx = 1:numel(prep_names)
    prep_idx = find(param.prep_names == prep_names(temp_prep_idx), 1);
    for epoch_idx = 1:numel(param.epoch_names)
        for cond_idx = 1:numel(param.cond_name_list)
            selected_idx = value_table.Preprocessing == param.prep_names(prep_idx) & value_table.EpochIndex == epoch_idx & value_table.ConditionIndex == cond_idx;
            selected_vals = value_table.(value_field)(selected_idx);
            selected_vals = selected_vals(isfinite(selected_vals));
            num_vals = numel(selected_vals);
            sd_vals = std(selected_vals, 'omitnan');
            quartiles = prctile(selected_vals, [25, 75]);

            temp_table = table(param.prep_names(prep_idx), param.prep_labels(prep_idx), param.prep_units(prep_idx), param.epoch_names(epoch_idx), epoch_idx, param.cond_name_list(cond_idx), cond_idx, ...
                num_vals, mean(selected_vals, 'omitnan'), sd_vals, sd_vals/sqrt(num_vals), median(selected_vals, 'omitnan'), quartiles(2)-quartiles(1), ...
                'VariableNames', {'Preprocessing', 'PreprocessingLabel', 'ValueUnit', 'Epoch', 'EpochIndex', 'Condition', 'ConditionIndex', 'N', 'Mean', 'SD', 'SEM', 'Median', 'IQR'});
            summary_table = [summary_table; temp_table];
        end
    end
end

end