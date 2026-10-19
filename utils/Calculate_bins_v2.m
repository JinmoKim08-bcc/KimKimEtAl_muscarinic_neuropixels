function [x_changed] = Calculate_bins_v2(x, fs, bin_time)

% x: column vector, n x 1
% fs: sampling rate of data x [Hz]
% bin_time: binning time [sec]

bin_pts = floor(bin_time * fs);
num_win = floor(numel(x) / bin_pts);

x_changed = ones(num_win, 1);

for win_idx = 1:num_win
    lb_idx = 1 + (win_idx-1) * bin_pts;
    ub_idx = win_idx * bin_pts;
    x_changed(win_idx) = sum(x(lb_idx:ub_idx))/bin_time;
end

end

