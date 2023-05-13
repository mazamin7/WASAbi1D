function [max_val, max_idx] = complex_max(varargin)
% This function returns the maximum and the index of the maximum
% between any number of complex scalar values.

% Initialize the maximum value and its index
max_val = varargin{1};
max_idx = 1;

% Loop through the input arguments
for i = 2:length(varargin)
    % Check if the current value is greater than the maximum value
    if abs(varargin{i}) > abs(max_val)
        % If it is, update the maximum value and its index
        max_val = varargin{i};
        max_idx = i;
    end
end
