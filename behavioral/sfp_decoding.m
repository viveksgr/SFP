%% Find raw decodability

nodor = 160;
wind = 7500; % Number of samples
dirs = {'C:\Work\SFP\sfp_behav_s01_correct';
    'C:\Work\SFP\sfp_behav_s02_correct';
    'C:\Work\SFP\sfp_behav_s04_correct'};
%   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));
corrmat_ = true;

svm_trainer2 = false;
svm_trainer3 = false;
svm_trainer4 = false;

%% Basic decoding
if corrmat_
corrmoda = zeros(3,1,2);
% corrmodb = zeros(3,1);
% corrmodp = zeros(3,1);
for ss = 1:length(dirs)
    load(fullfile(dirs{ss},'sfp_feats_main.mat'))
    Fless_mat = vertcat(fless_mat{:});
  % Fless_mat = vertcat(feat_mat{:});
     anatdir = fullfile('C:\Work\ARC\ARC\',sprintf('ARC%02d',ss),'single');

    
    if ss==3; s2 = 4; else; s2 = ss; end
    onsets = load(fullfile(anatdir,sprintf('conditions_NEMO%02d.mat',s2)),'onsets');
    onsets = onsets.onsets;
    group_vec = cell(nodor,1);
    unity = [];
    for ii2 = 1:nodor
        group_vec{ii2} = ii2*ones(length(onsets{ii2}),1);
        unity = blkdiag(unity,ones(length(onsets{ii2})));
    end
    group_vec = vertcat(group_vec{:});
    [~,argsort] = sort(vertcat(onsets{:}));
    group_vec = group_vec(argsort);
    unity = unity(argsort,argsort);
    utl_mask = logical(triu(ones(length(unity)),1)); % All possible odors
    
    Fless_mat_pruned = Fless_mat(:,1:wind);
    % Fless_mat_pruned = Fless_mat(:,[3 4 9:32]);
    Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
    Fless_corr = corrcoef(Fless_mat_pruned');
    

    % Trials belong to same odor
    M_on = logical(unity);
    M_on(logical(eye(size(M_on)))) = false; % Trivial correlation
   
    % Trials belong to different odors
    M_off = ~logical(unity);

    % Correlation across all trials belonging to different odors
    M_mid = Fless_corr;
    M_mid(logical(unity)) = nan;
    M_mid = nanmean(M_mid,2);
    
    % Correlation across all trials belonging to same odor
    M_mid_on = Fless_corr;
    M_mid_on(~M_on) = nan;
    M_mid_on = nanmean(M_mid_on,2);
   
    pval = ranksum( M_mid_on,M_mid)  ;  
    corrmoda(ss,1,1) = mean(Fless_corr(M_on));   
    corrmoda(ss,1,2) = mean(Fless_corr(M_off));   
    % corrmodb(ss) = std(M_mid_on-M_mid)/sqrt(4560);
    % corrmodp(ss)  = pval;
end
end

figure('Position',[0.5 0.5 320 240])
rsa_P1 = corrmoda;
S_mat = squeeze(mean(rsa_P1,1));
S_err = squeeze(std(rsa_P1,1))./sqrt(3)*1.96;
figure('Position',[0.5 0.5 400 250])
hold on
ngroups = size(S_mat, 1);
nbars = size(S_mat, 2);
bar(S_mat);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
x_m = [];
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/(2*nbars);
    errorbar(x, S_mat(:,i), S_err(:,i), 'k.');
    x_m = [x_m; x];
end

% legend({'Perceptual','Chemical','Mutual'})
% legend()
% Subject data points
c_s = {'r','g','b'}; % Data dots for subjects

for jj = 1:3
    plot([1 2],squeeze(rsa_P1(jj,1,:)),c_s{jj},'handle','off')
end
xticks([1 2])
xticklabels({'Same odor','Different odor'})
ylabel('Pattern correlation')


%% Decoding
nfolds = 10;
tic
numpcs = [13 11 11]; % 90% Variance

svm_trainer2 = true;
if svm_trainer2
corrmod = zeros(3,1);
predictions = cell(3,1);
baseline = cell(3,1);
pvalue = zeros(3,1);
figure()
corrmoda = zeros(3,1,2);
hold on
for ss = 1:length(dirs)
    load(fullfile(dirs{ss},'sfp_feats_main.mat')) 
    anatdir = fullfile('C:\Work\ARC\ARC\',sprintf('ARC%02d',ss),'single');
     
    if ss==3; s2 = 4; else; s2 = ss; end
    onsets = load(fullfile(anatdir,sprintf('conditions_NEMO%02d.mat',s2)),'onsets');
    onsets = onsets.onsets;
    group_vec = cell(nodor,1);
    unity = [];
    for ii2 = 1:nodor
        group_vec{ii2} = ii2*ones(length(onsets{ii2}),1);
        unity = blkdiag(unity,ones(length(onsets{ii2})));
    end
    group_vec = vertcat(group_vec{:});
    [~,argsort] = sort(vertcat(onsets{:}));
    group_vec = group_vec(argsort);
    unity = unity(argsort,argsort);
    utl_mask = logical(triu(ones(length(unity)),1)); % All possible odors
    
    % Fless_mat = vertcat(fless_mat{:});
    % Fless_mat_pruned = Fless_mat(:,1:100:wind);

    Fless_mat = vertcat(feat_mat{:});
    Fless_mat_pruned = Fless_mat(:,[3 4 9:31]);
    Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
    Fless_mat_pruned = zscore(Fless_mat_pruned,1);

    [coeff,Fless_mat_pruned,~,~,var] = pca(Fless_mat_pruned);
    Fless_mat_pruned = Fless_mat_pruned(:,1:numpcs(ss));
    
    subplot(1,3,ss)
    plot(cumsum(var))

    oid_ = 1:160;
    oid = oid_(group_vec)';
    
    [corrmod(ss),predictions_vec] = Classify_Permute_VS2(Fless_mat_pruned, oid, 5);

    % behavioral_corr
    accuracies = predictions_vec==oid;
    actual_behav = behav.behav(ss).ratings(group_vec,:);
    actual_behav = actual_behav(~accuracies ,:);
    predicted_behav = behav.behav(ss).ratings(predictions_vec,:);
    predicted_behav = predicted_behav(~accuracies ,:);
    predictions{ss} = iter_corr( actual_behav,predicted_behav);

    % Shuffle test on perceptual correlation on wrong trials
    baseline{ss}=  iter_corr_shuff( actual_behav,predicted_behav);
    corrmoda(ss,1,1) = mean( predictions{ss} );
    corrmoda(ss,1,2)= mean( baseline{ss} ,'all');

    T_shuff = mean(baseline{ss});
    t_stat = mean(predictions{ss});
    pValue(ss) = 2 * min(mean(T_shuff >= t_stat), mean(T_shuff <= t_stat));
end
end

figure('Position',[0 0 320 240]) 
hold on
bar(mean(corrmod))
errorbar(mean(corrmod),std(corrmod)./sqrt(3)*1.96)
hold on
c_s = {'r','g','b'};
for ss = 1:3; plot([1],corrmod(ss),c_s{ss},'Marker','.','MarkerSize',15); end 
yline(1/160)
ylabel('Performance')
savefig('svm')
print('svm','-dpng')
toc

p_value = arrayfun(@(x) ARC_computePValueOneTailed(x, 160, 4320),corrmod);
p_value_main = ARC_computePValueOneTailed(mean(corrmod), 160, 4320);


% Boxplots
figure()
VS_boxplot(predictions)
xticks(1:3)
xticklabels({'S1','S2','S3'})
ylabel('Performance (r)')
savefig('behav')
print('behav','-dpng')
p = [];
for ii = 1:3
     p(ii) = signrank(predictions{ss});
end
savefig('boxp')
print('boxp','-dpng')

clear unity fless_mat fless_mat_unn Fless_mat utl_mask
save('svmpred')
