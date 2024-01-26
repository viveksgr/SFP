function [dtw_values,p,fig] = SFP_computeDTW(T1, T2, optout,figs)
if nargin<4
    figs = false;
    fig = nan;
end
% Initialization
num_T1 = size(T1, 1); % number of observations in T1
num_T2 = size(T2, 1); % number of observations in T2
dtw_values = zeros(optout, 3); % Initialize the output matrix

% Create all possible pairs
all_pairs_T1 = nchoosek(1:num_T1, 2);
all_pairs_T2 = nchoosek(1:num_T2, 2);

% Randomly select pairs based on opt.out
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
    % dtw_values(i, 3) = dtw(T1(selected_pairs_T1(i, 1),:), T2(randi(num_T2),:));
    dtw_values(i, 3) = dtw(T1(selected_pairs_T1(i, 1),:), T2(selected_pairs_T2(i, 1),:));
end

[p1] = ranksum(dtw_values(:,1),dtw_values(:,3));
[p2] = ranksum(dtw_values(:,2),dtw_values(:,3));
p = max(p1,p2);

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
    boxplot(dtw_values)
    % bar([1 2 3],mean(dtw_values));
    hold on
    % errorbar([1 2 3],mean(dtw_values),1.96*std(dtw_values)./sqrt(optout),'.')
    xticks([1 2 3])
    xticklabels({'within trial (low)','within trial (high)', 'across trials'})
    xtickangle(45)
    ylim([min(dtw_values,[],"all"),max(dtw_values,[],"all")])
    ylabel('DTW distance')
    title(sprintf('Difference between traces pval: %.3f',p))
    savefig('SFP_dtw')
end
end
