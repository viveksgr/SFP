function [final_p_value, fig] = SFP_computeDTW_boot(T1, T2, optout, nboot, figs)

if nargin<5
    figs = false;
    fig = [];
end

% Initialization
num_T1 = size(T1, 1); % number of observations in T1
num_T2 = size(T2, 1); % number of observations in T2

% Initialize storage for bootstrap results
bootstrap_means = zeros(nboot, 3);

% Create all possible pairs
all_pairs_T1 = nchoosek(1:num_T1, 2);
all_pairs_T2 = nchoosek(1:num_T2, 2);

% Run bootstrap iterations
for b = 1:nboot
    
    dtw_values = zeros(optout, 3); % Initialize the output matrix
    
    % Randomly select pairs based on optout
    selected_pairs_T1 = all_pairs_T1(randperm(size(all_pairs_T1, 1), optout), :);
    selected_pairs_T2 = all_pairs_T2(randperm(size(all_pairs_T2, 1), optout), :);
    
    % Compute DTW for T1 Observations
    for i = 1:optout
        dtw_values(i, 1) = dtw(T1(selected_pairs_T1(i, 1),:), T1(selected_pairs_T1(i, 2),:));
    end

    % Compute DTW for T2 Observations
    for i = 1:optout
        dtw_values(i, 2) = dtw(T2(selected_pairs_T2(i, 1),:), T2(selected_pairs_T2(i, 2),:));
    end

    % Compute DTW Across T1 and T2 Observations
    for i = 1:optout
        dtw_values(i, 3) = dtw(T1(selected_pairs_T1(i, 1),:), T2(selected_pairs_T2(i, 1),:));
    end
    
    % Store the means in the bootstrap storage
    bootstrap_means(b, :) = mean(dtw_values);
end

p_value_T1 = bootp_val(bootstrap_means(:, 3),bootstrap_means(:, 1));
p_value_T2 = bootp_val(bootstrap_means(:, 3),bootstrap_means(:, 2));

% Take the maximum p-value as the final p-value
final_p_value = max([p_value_T1, p_value_T2]);


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
    boxplot(bootstrap_means)
    % bar([1 2 3],mean(dtw_values));
    hold on
    % errorbar([1 2 3],mean(dtw_values),1.96*std(dtw_values)./sqrt(optout),'.')
    xticks([1 2 3])
    xticklabels({'within trial (low)','within trial (high)', 'across trials'})
    xtickangle(45)
    ylim([min(bootstrap_means,[],"all"),max(bootstrap_means,[],"all")])
    ylabel('DTW distance')
    title(sprintf('Difference between traces pval: %.3f',final_p_value))
    savefig('SFP_dtw')
end


end

% Compute p_val
function p_val = bootp_val(t1,t2)

% Compute the observed test statistic t_sq
t_sq = (mean(t1) - mean(t2)) / sqrt(var(t1) + var(t2));

% Create null distribution by shifting the means
grp1 = t1 - mean(t1) + mean([t1; t2]); % For group 1
grp2 = t2 - mean(t2) + mean([t1; t2]); % For group 2

% Compute the test statistic for each bootstrap sample to form the null distribution
tdist = (grp1 - grp2) ./ sqrt(var(grp1) + var(grp2));

% Compute the two-tailed p-value
p_val = 2 * min(sum(tdist >= t_sq) / length(tdist), sum(tdist <= t_sq) / length(tdist));
end
