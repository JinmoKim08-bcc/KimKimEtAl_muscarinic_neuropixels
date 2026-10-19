function fig = Plot_distribution_figure(value_table, summary_table, pairwise_table, prep_name, value_field, variable_label, is_unsigned, do_crop, param, panel_rows)
% panel_rows: rows to draw (1: mean +/- SEM bar, 2: histogram, 3: cumulative histogram)

num_epochs = numel(param.epoch_names);
num_cond = numel(param.cond_name_list);
prep_idx = find(param.prep_names == prep_name, 1);
value_unit = param.prep_units(prep_idx);
hist_bin_width = param.hist_bin_widths(prep_idx);
display_order = param.bar_cond_order(:);

% |RR-LL| of TrialAvgZScore in [0, 2]: one patch per bar and fixed y-axis
is_main_fig = do_crop && prep_name == "TrialAvgZScore" && value_field == "AbsRRMinusLLValue";

fig = figure('Color', 'w', 'Position', [30, 30, 390*num_epochs, 300*numel(panel_rows)]);
tile_layout = tiledlayout(fig, numel(panel_rows), num_epochs, 'TileSpacing', 'compact', 'Padding', 'compact');

for epoch_idx = 1:num_epochs
    base_idx = value_table.Preprocessing == prep_name & value_table.EpochIndex == epoch_idx;

    %%% Mean +/- SEM bar plot
    if ismember(1, panel_rows)
        ax = nexttile(tile_layout, (find(panel_rows == 1, 1)-1)*num_epochs+epoch_idx);
        hold(ax, 'on');
        selected_summary = summary_table(summary_table.Preprocessing == prep_name & summary_table.EpochIndex == epoch_idx, :);
        [~, row_order] = ismember(display_order, selected_summary.ConditionIndex);
        selected_summary = selected_summary(row_order, :);

        if is_main_fig
            % separate patches so that each bar can be recolored in Illustrator
            bar_width = 0.75;
            for bar_idx = 1:num_cond
                bar_mean = selected_summary.Mean(bar_idx);
                patch(ax, [bar_idx-bar_width/2, bar_idx+bar_width/2, bar_idx+bar_width/2, bar_idx-bar_width/2], [0, 0, bar_mean, bar_mean], param.cond_colors(display_order(bar_idx), :), 'EdgeColor', 'none', 'FaceAlpha', 1, 'HandleVisibility', 'off');
            end
        else
            bar_handle = bar(ax, 1:num_cond, selected_summary.Mean, 0.75, 'FaceColor', 'flat', 'EdgeColor', 'none');
            bar_handle.CData = param.cond_colors(display_order, :);
        end
        errorbar(ax, 1:num_cond, selected_summary.Mean, selected_summary.SEM, 'k', 'LineStyle', 'none', 'LineWidth', 1.2, 'CapSize', 8);
        if ~is_unsigned
            yline(ax, 0, 'k--', 'HandleVisibility', 'off');
        end
        title(ax, param.epoch_names(epoch_idx), 'Interpreter', 'none');
        xlim(ax, [0.4, num_cond+0.6]);
        xticks(ax, 1:num_cond);
        xticklabels(ax, param.cond_name_list(display_order));
        if is_unsigned
            ylabel(ax, 'Mean difference in normalized firing rate');
        else
            ylabel(ax, 'Mean normalized firing rate');
        end

        % brackets with Bonferroni-adjusted p-values (short pairs first)
        selected_pairwise = pairwise_table(pairwise_table.Preprocessing == prep_name & pairwise_table.EpochIndex == epoch_idx, :);
        bracket_pairs = nchoosek(1:num_cond, 2);
        [~, sort_idx] = sort(diff(bracket_pairs, 1, 2));
        bracket_pairs = bracket_pairs(sort_idx, :);
        y_top = max([0; selected_summary.Mean+selected_summary.SEM], [], 'omitnan');
        y_bottom = min([0; selected_summary.Mean-selected_summary.SEM], [], 'omitnan');
        y_span = y_top - y_bottom;
        bracket_base = y_top + 0.08*y_span;
        bracket_gap = 0.17*y_span;
        bracket_cap = 0.045*y_span;

        for pair_idx = 1:size(bracket_pairs, 1)
            x_pair = bracket_pairs(pair_idx, :);
            name1 = param.cond_name_list(display_order(x_pair(1)));
            name2 = param.cond_name_list(display_order(x_pair(2)));
            selected_idx = (selected_pairwise.Group1 == name1 & selected_pairwise.Group2 == name2) | (selected_pairwise.Group1 == name2 & selected_pairwise.Group2 == name1);
            p_adj = selected_pairwise.PBonferroniWithinEpoch(find(selected_idx, 1));
            if p_adj < 1e-4
                p_str = 'p=<10^{-4}';
            else
                p_str = sprintf('p=%.3g', p_adj);
            end

            bracket_y = bracket_base + (pair_idx-1)*bracket_gap;
            plot(ax, [x_pair(1), x_pair(1), x_pair(2), x_pair(2)], [bracket_y, bracket_y+bracket_cap, bracket_y+bracket_cap, bracket_y], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
            text(ax, mean(x_pair), bracket_y+1.15*bracket_cap, p_str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', 'k');
        end
        ylim(ax, [y_bottom-0.04*y_span, bracket_base+(size(bracket_pairs, 1)-1)*bracket_gap+bracket_cap+0.16*y_span]);

        if is_main_fig
            bar_ymax = [0.4, 0.2, 0.3, 0.5];
            bar_ystep = [0.1, 0.05, 0.05, 0.05];
            ylim(ax, [0, bar_ymax(epoch_idx)]);
            yticks(ax, 0:bar_ystep(epoch_idx):bar_ymax(epoch_idx));
        end
        grid(ax, 'on');
        set(ax, 'FontSize', 10, 'LineWidth', 1, 'Box', 'off', 'TickDir', 'out');
    end

    %%% Histogram bin edges
    if do_crop
        edges = linspace(param.abs_diff_xlim(1), param.abs_diff_xlim(2), round(diff(param.abs_diff_xlim)/hist_bin_width)+1);
    else
        selected_vals = value_table.(value_field)(base_idx);
        selected_vals = selected_vals(isfinite(selected_vals));
        edge_min = floor(min(selected_vals)/hist_bin_width)*hist_bin_width;
        edge_max = ceil(max(selected_vals)/hist_bin_width)*hist_bin_width;
        edges = edge_min:hist_bin_width:edge_max;
        if edges(end) <= max(selected_vals)
            edges(end+1) = edges(end) + hist_bin_width;
        end
    end

    %%% Probability histogram
    if ismember(2, panel_rows)
        ax = nexttile(tile_layout, (find(panel_rows == 2, 1)-1)*num_epochs+epoch_idx);
        hold(ax, 'on');
        hist_handles = gobjects(num_cond, 1);
        for cond_idx = 1:num_cond
            selected_vals = value_table.(value_field)(base_idx & value_table.ConditionIndex == cond_idx);
            selected_vals = selected_vals(isfinite(selected_vals));
            hist_handles(cond_idx) = histogram(ax, selected_vals, 'BinEdges', edges, 'Normalization', 'probability', 'DisplayStyle', param.hist_display_style, ...
                'FaceColor', param.cond_colors(cond_idx, :), 'FaceAlpha', param.hist_face_alpha, 'EdgeColor', param.cond_colors(cond_idx, :), 'EdgeAlpha', param.hist_edge_alpha, 'LineWidth', 0.35);
        end
        if ~is_unsigned
            xline(ax, 0, 'k--', 'HandleVisibility', 'off');
        end
        xlim(ax, [edges(1), edges(end)]);
        xlabel(ax, sprintf('%s [%s]', variable_label, value_unit));
        ylabel(ax, 'Fraction of neurons');
        if do_crop
            % number of units above the x range for each condition
            num_excluded = zeros(1, num_cond);
            excluded_text = strings(1, num_cond);
            for bar_idx = 1:num_cond
                selected_vals = value_table.(value_field)(base_idx & value_table.ConditionIndex == display_order(bar_idx));
                num_excluded(bar_idx) = sum(isfinite(selected_vals) & selected_vals > param.abs_diff_xlim(2));
                temp_name = char(param.cond_name_list(display_order(bar_idx)));
                excluded_text(bar_idx) = sprintf('%s %d', temp_name(1), num_excluded(bar_idx));
            end
            title(ax, {sprintf('Histogram (bin = %.2f %s)', hist_bin_width, value_unit); sprintf('Not shown (>%.2g): %s | total %d', param.abs_diff_xlim(2), char(strjoin(excluded_text, ' | ')), sum(num_excluded))}, 'Interpreter', 'none');
        else
            title(ax, sprintf('Histogram (bin = %.2f %s)', hist_bin_width, value_unit));
        end
        grid(ax, 'on');
        set(ax, 'FontSize', 10, 'LineWidth', 1, 'Box', 'off', 'TickDir', 'out');
        if epoch_idx == 1
            legend(ax, hist_handles, cellstr(param.cond_name_list), 'Location', 'best', 'Box', 'off');
        end
    end

    %%% Cumulative histogram from the same bin counts
    if ismember(3, panel_rows)
        ax = nexttile(tile_layout, (find(panel_rows == 3, 1)-1)*num_epochs+epoch_idx);
        hold(ax, 'on');
        for cond_idx = 1:num_cond
            selected_vals = value_table.(value_field)(base_idx & value_table.ConditionIndex == cond_idx);
            selected_vals = selected_vals(isfinite(selected_vals));
            hist_counts = histcounts(selected_vals, edges);
            if do_crop
                % units above the x range stay in the denominator
                cum_fraction = cumsum(hist_counts)/numel(selected_vals);
            else
                cum_fraction = cumsum(hist_counts)/sum(hist_counts);
            end
            stairs(ax, edges, [0, cum_fraction], 'Color', param.cond_colors(cond_idx, :), 'LineWidth', 1.8);
        end
        if ~is_unsigned
            xline(ax, 0, 'k--', 'HandleVisibility', 'off');
        end
        xlim(ax, [edges(1), edges(end)]);
        ylim(ax, [0, 1]);
        xlabel(ax, sprintf('%s [%s]', variable_label, value_unit));
        ylabel(ax, 'Cumulative fraction');
        title(ax, 'Cumulative histogram');
        grid(ax, 'on');
        set(ax, 'FontSize', 10, 'LineWidth', 1, 'Box', 'off', 'TickDir', 'out');
    end
end

if do_crop
    range_title = sprintf('x range [%.2g, %.2g]; overflow not shown', param.abs_diff_xlim(1), param.abs_diff_xlim(2));
else
    range_title = 'full x range';
end
sgtitle(tile_layout, sprintf('%s | %s | pooled neurons | %s', variable_label, param.prep_labels(prep_idx), range_title), 'Interpreter', 'none');

end
