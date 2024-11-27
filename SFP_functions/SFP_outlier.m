function T1_clean = SFP_outlier(T1,threshold)
% Let's say T1 is your data matrix for the first group, of size [num_trials x num_timepoints]

if nargin<2
    threshold = 3;
end

% Compute mean and std across trials
mean_T1 = mean(T1, 1);
std_T1 = std(T1, 0, 1);

% Compute Z-scores
Z_T1 = abs((T1 - mean_T1) ./ std_T1);

% Identify outliers
outliers = any(Z_T1 > threshold, 2); % threshold might be something like 2 or 3

% Remove outliers
T1_clean = T1(~outliers, :);
