function [gridFPR, TPRmat] = Interpolate_ROC_curve_v0b(fpr_list, tpr_list, param)

N = numel(fpr_list);
    allFPR  = cat(1, fpr_list{:});
    gridFPR = unique(allFPR);                 

    TPRmat = nan(N, numel(gridFPR));
    for i = 1:N
        fpr = fpr_list{i}(:);
        tpr = tpr_list{i}(:);

        [fpr_sorted, idx] = sort(fpr);
        tpr_sorted        = tpr(idx);
        [fpr_u, ia] = unique(fpr_sorted, 'stable');
        tpr_u       = tpr_sorted(ia);

        TPRmat(i,:) = interp1(fpr_u, tpr_u, gridFPR, 'linear', 0);
    end
end