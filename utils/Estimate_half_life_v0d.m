function [t_half, R2] = Estimate_half_life_v0d(time, y, min_use_y)

valid = (y > min_use_y);

time = time(valid);
y = y(valid);

ln_y = log(y);

% assumption: y follows A*exp(-t/t_half)
% linear regression
p = polyfit(time, ln_y, 1);

% p(1): slope, p(2): offset
slope = p(1);
% ln y = ln A + (–λ) t => slope = –λ
lambda_est = -slope;

% tau decay
t_half = log(2) / lambda_est;

% estimate
ln_y_pred = polyval(p, time);
y_pred    = exp(ln_y_pred)';

% R^2 calculation
SS_res = sum((y - y_pred).^2);
SS_tot = sum((y - mean(y)).^2);
R2      = 1 - SS_res / SS_tot;

end

