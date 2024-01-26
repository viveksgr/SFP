function [p_val, actual_statistic, permuted_stats,fig] = SFP_computeDTW_perm(T1, T2, optout, n_perms, figs)

if nargin<5
    figs = false;
    fig = [];
end

% Combine the two groups T1 and T2
combined_T = [T1; T2];
num_T1 = size(T1, 1);
num_T2 = size(T2, 1);
total_samples = num_T1 + num_T2;

% % Compute the actual mean DTW distance between T1 and T2
% actual_statistic_mat = zeros(n_perms, 1);
% for i = 1:n_perms  
%     % Compute the permuted mean DTW distance
%     actual_statistic_mat(i) = mean(SFP_computeDTW_helper(T1,T2, optout));
% end
% actual_statistic = median(actual_statistic_mat);
actual_statistic = mean(computeDTW(T1,T2, optout));

% Initialize array to store permutation statistics
permuted_stats = zeros(n_perms, 1);

% Run permutation tests
for i = 1:n_perms
    % Shuffle the labels
    shuffled_indices = randperm(total_samples);
    new_T1 = combined_T(shuffled_indices(1:num_T1), :);
    new_T2 = combined_T(shuffled_indices(num_T1+1:end), :);
    
    % Compute the permuted mean DTW distance
    permuted_stats(i) = mean(computeDTW(new_T1, new_T2, optout));
end

% Compute the p-value
p_val = sum(abs(permuted_stats) >= abs(actual_statistic)) / n_perms;


if figs
    fig = figure('Position',[0.1 0.1 1280 480]);
    nplots = min([50,size(T1,1),size(T2,1)]);
    hold on
    subplot(1,2,1)
    hold on
    plot(mean(T1)','Linewidth',2,'color',[0.8500 0.3250 0.0980])
    plot(mean(T2)','Linewidth',2,'color',[0 0.4470 0.741])
    legend({'T1','T2'})
    plot(T1(1:nplots ,:)','color',[0.8500 0.3250 0.0980 0.2],'HandleVisibility','off')
    plot(T2(1:nplots ,:)','color',[0 0.4470 0.741 0.2],'HandleVisibility','off')
    ylabel('Amplitude')
    xlabel('Time')


    subplot(1,2,2)
    hold on
    histogram(permuted_stats)
    xline(actual_statistic)
    % bar([1 2 3],mean(dtw_values));
    hold on
    % errorbar([1 2 3],mean(dtw_values),1.96*std(dtw_values)./sqrt(optout),'.')
    % xticks([1 2 3])
    % xticklabels({'within trial (low)','within trial (high)', 'across trials'})
    % xtickangle(45)
    % ylim([min(permuted_stats,[],"all"),max(permuted_stats,[],"all")])
    ylabel('DTW distance histogram')
    title(sprintf('Difference between traces pval: %.3f',p_val))
    savefig('SFP_dtw')
end



end
% Compute the two-tailed p-value
% p_val_two_tailed = sum(abs(permuted_stats) >= abs(actual_statistic) | abs(permuted_stats) <= -abs(actual_statistic)) / n_perms;

% Helper function to compute DTW
function [mean_dtw_distance] = computeDTW(T1, T2, optout)
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
