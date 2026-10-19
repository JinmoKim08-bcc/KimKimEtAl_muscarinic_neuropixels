function x_smooth = Perform_smoothing_fr(x, param)

x_smooth = zeros(size(x));

for unit_idx = 1:size(x, 1)
    temp_x = conv(x(unit_idx, :), param.smoothing_window, 'same');

    x_smooth(unit_idx, :) = temp_x;

end


end

