function [ val ] = weightedAverage(weights, values)
%WEIGHTEDAVERAGE Calculate the weighted average of 'values' by applying 
% the 'weights'
%
%   values - Data points to average, one per row.
%  weights - Weight to apply to each data point, one per row.
%
%  Returns:
%     val  - The weighted average of 'values'.

    % Apply the weights to the values by taking the dot-product between the
    % two vectors.
    val = weights' * values;
    
    % Divide the sum by the weights
    val = val ./ sum(weights, 1);
end

