function [sig_idx, latency] = Find_response_latency(x_baseline, x, max_pts)
    % input args
    % x_baseline: num_units x time_points (baseline activity)
    % x: num_units x time_points
    % max_pts: considered maximal time points for the response

    % Find significant latency based on the paper, Inagaki, et al., Cell, 2022.
    x_baseline_norm = x_baseline;
    x_baseline_norm = (x_baseline_norm - mean(x_baseline_norm, 2)) ./ std(x_baseline_norm, [], 2);
    x_baseline_norm(x_baseline_norm<0) = 0;
    
    num_units = size(x_baseline, 1);

    %%% tp0001
    alpha = 0.001;
    [lambdahat, lambdaci] = poissfit(x_baseline_norm', alpha);
    thresh0001 = lambdaci(2, :)';
    tp0001 = zeros(num_units, 1);
    for unit_idx = 1:num_units
        if numel(find(x(unit_idx, :) > thresh0001(unit_idx), 1, 'first')) ~= 0
            tp0001(unit_idx) = find(x(unit_idx, :) > thresh0001(unit_idx), 1, 'first');
        else
            tp0001(unit_idx) = 0;
        end
    end
    
    %%% tp005
    alpha = 0.05;
    [lambdahat, lambdaci] = poissfit(x_baseline_norm', alpha);
    thresh005 = lambdaci(2, :)';

    %%% Find latency
    latency = zeros(num_units, 1);
    for unit_idx = 1:num_units
        if numel(find(x(unit_idx, :) > thresh0001(unit_idx), 1, 'first')) ~= 0
            latency(unit_idx) = find(x(unit_idx, 1:tp0001(unit_idx)) > thresh005(unit_idx), 1, 'last');
        else
            latency(unit_idx) = 0;
        end
    end
    [M, I] = sort(latency, 'ascend');
    sig_idx = I(M~=0 & M<=max_pts);
    latency = latency(I(M~=0 & M<=max_pts));
    
end

