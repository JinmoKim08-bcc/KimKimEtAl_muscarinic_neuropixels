function [results_corr] = Calculate_population_pearson_corr_fast(x1, x2)

% x1: N x T
% x2: N x T

[N, T] = size(x1);

% x1, x2: N x T
R = corrcoef([x1, x2]);  % (T+T)×(T+T) correlation matrix
results_corr = R(1:T, T+1:2*T);

end

