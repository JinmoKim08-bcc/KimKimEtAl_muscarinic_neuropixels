function [omnibus_table, pairwise_table] = Calculate_pooled_statistics(value_table, param, value_field, prep_names)

omnibus_table = table();
pairwise_table = table();
num_pairs = size(param.cond_pairs, 1);

for temp_prep_idx = 1:numel(prep_names)
    prep_idx = find(param.prep_names == prep_names(temp_prep_idx), 1);
    for epoch_idx = 1:numel(param.epoch_names)
        cond_vals = cell(numel(param.cond_name_list), 1);
        all_vals = [];
        all_groups = [];
        for cond_idx = 1:numel(param.cond_name_list)
            selected_idx = value_table.Preprocessing == param.prep_names(prep_idx) & value_table.EpochIndex == epoch_idx & value_table.ConditionIndex == cond_idx;
            selected_vals = value_table.(value_field)(selected_idx);
            cond_vals{cond_idx} = selected_vals(isfinite(selected_vals));
            all_vals = [all_vals; cond_vals{cond_idx}];
            all_groups = [all_groups; repmat(cond_idx, numel(cond_vals{cond_idx}), 1)];
        end

        %%% Kruskal-Wallis test across conditions
        p_kw = kruskalwallis(all_vals, all_groups, 'off');

        temp_table = table(param.prep_names(prep_idx), param.prep_labels(prep_idx), param.prep_units(prep_idx), param.epoch_names(epoch_idx), epoch_idx, p_kw, NaN, ...
            'VariableNames', {'Preprocessing', 'PreprocessingLabel', 'ValueUnit', 'Epoch', 'EpochIndex', 'PRaw', 'PHolmAcrossEpochs'});
        omnibus_table = [omnibus_table; temp_table];

        %%% Wilcoxon rank-sum test and Cliff's delta for each pair
        for pair_idx = 1:num_pairs
            cond1_idx = param.cond_pairs(pair_idx, 1);
            cond2_idx = param.cond_pairs(pair_idx, 2);
            vals1 = cond_vals{cond1_idx};
            vals2 = cond_vals{cond2_idx};

            p_raw = ranksum(vals1, vals2);

            ranks = tiedrank([vals1; vals2]);
            U = sum(ranks(1:numel(vals1))) - numel(vals1)*(numel(vals1)+1)/2;
            cliffs_delta = 2*U/(numel(vals1)*numel(vals2)) - 1;

            median1 = median(vals1, 'omitnan');
            median2 = median(vals2, 'omitnan');

            temp_table = table(param.prep_names(prep_idx), param.prep_labels(prep_idx), param.prep_units(prep_idx), param.epoch_names(epoch_idx), epoch_idx, ...
                param.cond_name_list(cond1_idx), param.cond_name_list(cond2_idx), numel(vals1), numel(vals2), median1, median2, median1-median2, cliffs_delta, p_raw, NaN, ...
                'VariableNames', {'Preprocessing', 'PreprocessingLabel', 'ValueUnit', 'Epoch', 'EpochIndex', 'Group1', 'Group2', 'N1', 'N2', 'Median1', ...
                'Median2', 'MedianDifference', 'CliffsDelta', 'PRaw', 'PBonferroniWithinEpoch'});
            pairwise_table = [pairwise_table; temp_table];
        end
    end
end

%%% Holm correction across epochs, Bonferroni correction across pairs within each epoch
for temp_prep_idx = 1:numel(prep_names)
    selected_idx = omnibus_table.Preprocessing == prep_names(temp_prep_idx);
    [p_sorted, sort_idx] = sort(omnibus_table.PRaw(selected_idx), 'ascend');
    num_tests = numel(p_sorted);
    p_holm = zeros(num_tests, 1);
    p_holm(sort_idx) = min(cummax((num_tests-(1:num_tests)'+1).*p_sorted), 1);
    omnibus_table.PHolmAcrossEpochs(selected_idx) = p_holm;

    for epoch_idx = 1:numel(param.epoch_names)
        selected_idx = pairwise_table.Preprocessing == prep_names(temp_prep_idx) & pairwise_table.EpochIndex == epoch_idx;
        pairwise_table.PBonferroniWithinEpoch(selected_idx) = min(num_pairs*pairwise_table.PRaw(selected_idx), 1);
    end
end

end