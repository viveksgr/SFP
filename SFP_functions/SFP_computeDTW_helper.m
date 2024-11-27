
% Helper function to compute DTW
function [mean_dtw_distance] = SFP_computeDTW_helper(T1, T2, optout)
    num_T1 = size(T1, 1); % number of observations in T1
    num_T2 = size(T2, 1); % number of observations in T2
    dtw_distances = zeros(optout, 1); % Initialize the output array
    
    % Randomly select pairs for comparison
    for i = 1:optout
        dtw_distances(i) = dtw(T1(randi(num_T1), :), T2(randi(num_T2), :));
    end
    
    % Compute the mean DTW distance
    mean_dtw_distance = mean(dtw_distances);
end
