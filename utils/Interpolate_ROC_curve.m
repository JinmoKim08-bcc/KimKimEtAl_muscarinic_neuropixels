function [fpr, tpr] = Interpolate_ROC_curve(fpr_list, tpr_list, param)

allFPR = cat(1, fpr_list{:});
gridFPR = unique(allFPR);

allTPR = nan(size(fpr_list, 2), numel(gridFPR));

for i = 1:size(fpr_list, 2)
    fpr = fpr_list{i}(:);
    tpr = tpr_list{i}(:);

    [fpr_sorted, sortIdx] = sort(fpr);
    tpr_sorted = tpr(sortIdx);

    [fpr_u, ia] = unique(fpr_sorted, 'stable');
    tpr_u       = tpr_sorted(ia);

    allTPR(i,:) = interp1(fpr_u, tpr_u, gridFPR, 'linear', 0);
end

fpr = gridFPR;
tpr = allTPR;

end

