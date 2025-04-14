%% Basic decoding

numpcs = [10 10 10]; % 90% Variance

corrmod = zeros(3,1);
predictions = cell(3,1);
baseline = cell(3,1);
pvalue = zeros(3,1);
figure()
corrmoda = zeros(3,1,2);
hold on

for ss = 1:4

    Fless_mat_pruned  = feat_mat{ss}(:,[3 4 9:21 23:31]);
    % Fless_mat_pruned = Fless_mat(:,[3 4 9:31]);

    Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
    Fless_mat_pruned = zscore(Fless_mat_pruned,1);

    [coeff,Fless_mat_pruned,~,~,var] = pca(Fless_mat_pruned);
    cumvar = cumsum(var);
    numpc = sum(cumvar<70)+1
    Fless_mat_pruned = Fless_mat_pruned(:,1: numpc);
    
    oid = odor_id{ss};

    % % oid = oid(randperm(length(oid)));
    % if ss==4
    %     oid = circshift(oid,1);
    % end
    
    [~,predictions_vec] = Classify_Permute_VS2(Fless_mat_pruned, oid, 200);
      % [~,predictions_vec] = SFP_regress_nested2_normed(Fless_mat_pruned, oid, 20);
      % [~,predictions_vec] =  ClassifyNestedPCA_SVM(Fless_mat_pruned, oid, 100);

    % behavioral_corr
    accuracies = predictions_vec==oid;
    
    % Accuracies
    corrmod(ss) = sum(accuracies)/length(accuracies);
end


p_value_svm = arrayfun(@(x) ARC_computePValueOneTailed(x, 10, 200)/2,corrmod)

% Basic SVM
figure('Position',[0 0 320 240]) 
hold on
bar(mean(corrmod))
errorbar(mean(corrmod),std(corrmod)./sqrt(3)*1.96)
hold on
c_s = {'r','g','b','m'};
for ss = 1:4; plot([1],corrmod(ss),c_s{ss},'Marker','.','MarkerSize',15); end 
yline(1/10)
ylabel('Performance')
xticks(1)
xticklabels({''})
% savefig(fullfile(savepath,'svm'))
% print(fullfile(savepath,'svm'),'-dpng')

% p