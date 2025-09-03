%% Basic decoding
% Load 

nS = 3;
corrmod = zeros(nS,1);
predictions = cell(nS,1);
baseline = cell(nS,1);
pvalue = zeros(nS,1);
figure()
corrmoda = zeros(nS,1,2);
hold on

for ss = 1:3
    fprintf('Subject %02d\n',ss)

    Fless_mat_pruned  = feat_mat{ss}(:,[3 4 9:21 23:31]);
    % Fless_mat_pruned = Fless_mat(:,[3 4 9:31]);

    Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
    Fless_mat_pruned = (zscore(Fless_mat_pruned,1));
  
     % Fless_mat_pruned = vs_normalizer(Fless_mat_pruned);

    [coeff,Fless_mat_pruned,~,~,var] = pca(Fless_mat_pruned);
    cumvar = cumsum(var);
    numpc = sum(cumvar<70)+1;
    Fless_mat_pruned = Fless_mat_pruned(:,1: numpc);
    
    oid = odor_id{ss};

    % % oid = oid(randperm(length(oid)));
    % if ss==4
    %     oid = circshift(oid,1);
    % end
    
    [~,predictions_vec] = Classify_Permute_VS2(Fless_mat_pruned, oid, 200);
      % [~,predictions_vec] = SFP_regress_nested2_normed(Fless_mat_pruned, oid, 20);
      % [~,predictions_vec] =  ClassifyNestedPCA_SVM(Fless_mat_pruned, oid, 100);
    % [~,predictions_vec] = SFP_regress_nested2_normed_pca(Fless_mat_pruned, oid, 20);

    % behavioral_corr
    accuracies = predictions_vec==oid;
    
    % Accuracies
    corrmod(ss) = sum(accuracies)/length(accuracies);
end


p_value_svm = arrayfun(@(x) ARC_computePValueOneTailed(x, 10, 200),corrmod)

% Basic SVM
figure('Position',[0 0 320 240]) 
hold on
bar(mean(corrmod))
errorbar(mean(corrmod),std(corrmod)./sqrt(3)*1.96)
hold on
c_s = {'r','g','b','m'};
for ss = 1:nS; plot([1],corrmod(ss),c_s{ss},'Marker','.','MarkerSize',15); end 
yline(1/10)
ylabel('Performance')
xticks(1)
xticklabels({''})
% savefig(fullfile(savepath,'svm'))
% print(fullfile(savepath,'svm'),'-dpng')

% %% Same vs different 
% wind = 7500;
% corrmoda = zeros(4,1,2); % Subject-wise pattern correlations
% sig_clean = true;
% for ss = 1:4
% 
%     Fless_mat = fless_mat{ss};
%     % Fless_mat = vertcat(feat_mat{:});
% 
%     unity = 1-odor_id{1}+odor_id{1}';
%     unity(unity~=1) = 0;
% 
%     utl_mask = logical(triu(ones(length(unity)),1)); % All possible odors
%     Fless_mat_pruned = Fless_mat(:,1:wind);
% 
%     % Fless_mat_pruned = Fless_mat(:,[3 4 9:32]);
%     Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
%     if sig_clean
%         Fless_mat_pruned = Fless_mat_pruned-mean(Fless_mat_pruned);
%     end
%     Fless_corr = corrcoef(Fless_mat_pruned');
% 
%     % Trials belong to same odor
%     M_on = logical(unity);
%     M_on(logical(eye(size(M_on)))) = false; % Trivial correlation
% 
%     % Trials belong to different odors
%     M_off = ~logical(unity);
% 
%     % Correlation across all trials belonging to different odors
%     M_mid = Fless_corr;
%     M_mid(logical(unity)) = nan;
%     M_mid = nanmean(M_mid,2);
% 
%     % Correlation across all trials belonging to same odor
%     M_mid_on = Fless_corr;
%     M_mid_on(~M_on) = nan;
%     M_mid_on = nanmean(M_mid_on,2);
% 
%     pval = ranksum( M_mid_on,M_mid)  ;
%     corrmoda(ss,1,1) = mean(Fless_corr(M_on));
%     corrmoda(ss,1,2) = mean(Fless_corr(M_off));
%     % corrmodb(ss) = std(M_mid_on-M_mid)/sqrt(4560);
%     corrmodp(ss)  = pval;
% end
% 
% % Figure plotting
% figure('Position',[0.5 0.5 320 240])
% rsa_P1 = corrmoda;
% S_mat = squeeze(mean(rsa_P1,1));
% S_err = squeeze(std(rsa_P1,1))./sqrt(3)*1.96;
% figure('Position',[0.5 0.5 400 250])
% hold on
% ngroups = size(S_mat, 1);
% nbars = size(S_mat, 2);
% bar(S_mat);
% % Calculating the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% x_m = [];
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/(2*nbars);
%     errorbar(x, S_mat(:,i), S_err(:,i), 'k.');
%     x_m = [x_m; x];
% end
% c_s = {'r','g','b','m'}; % Data dots for subjects
% for jj = 1:4
%     plot([1 2],squeeze(rsa_P1(jj,1,:)),c_s{jj},'handle','off')
% end
% xticks([1 2])
% xticklabels({'Same odor','Different odor'})
% ylabel('Pattern correlation')
