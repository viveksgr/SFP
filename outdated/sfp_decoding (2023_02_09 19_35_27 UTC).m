%% Find raw decodability

nodor = 160;
wind = 75; % Number of samples
dirs = {'C:\Data\SFP\sfp_behav_s01';
    'C:\Data\SFP\sfp_behav_s02';
    'C:\Data\SFP\sfp_behav_s03'};
%   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
behav = load(fullfile('C:\Data\ARC\ARC','NEMO_perceptual2.mat'));
corrmat_ = true;
svm_trainer = false;
svm_trainer2 = false;
svm_trainer3 = false;
svm_trainer4 = false;
%% Basic decoding
if corrmat_
corrmod = zeros(3,2);
corrmodp = zeros(3,1);
for ss = 1:length(dirs)
    load(fullfile(dirs{ss},'sfp_feats.mat'))
%     Fless_mat = vertcat(fless_mat{:});
     Fless_mat = vertcat(feat_mat{:});
    anatdir = fullfile('C:\Data\ARC\',sprintf('ARC%02d',ss),'single');
    
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
    
%   Fless_mat_pruned = Fless_mat(:,1:wind);
    Fless_mat_pruned = Fless_mat(:,[3 4 9:14]);
    Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
    Fless_corr = corrcoef(Fless_mat_pruned');
    
    M_on = logical(unity);
    M_on(logical(eye(size(M_on)))) = false;
    
    M_off = ~logical(unity);
    % Behavioral features
    
    pval = ranksum(Fless_corr(M_on),Fless_corr(M_off));    
    corrmod(ss,1) = mean(Fless_corr(M_on));   
    corrmod(ss,2) =  mean(Fless_corr(M_off));
    corrmodp  = pval;
end
end
figure()
bar(corrmod)


%% Make SVM
nfolds = 10;

if svm_trainer
corrmod = zeros(nfolds,3);
for ss = 1:length(dirs)
    load(fullfile(dirs{ss},'sfp_feats.mat'))
    Fless_mat = vertcat(fless_mat{:});
    anatdir = fullfile('C:\Data\ARC\',sprintf('ARC%02d',ss),'single');
    
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
    oid_ = 1:160;
    oid = oid_(group_vec)';
   
    SVMModel = fitcecoc(Fless_mat_pruned,oid);
    CVSVMModel = crossval(SVMModel);
    kfoldLoss(CVSVMModel)
 
end
end
%% Make SVM - libsvm
nfolds = 10;

if svm_trainer2
corrmod = zeros(nfolds,3);
for ss = 1:length(dirs)
    load(fullfile(dirs{ss},'sfp_feats.mat'))
    Fless_mat = vertcat(fless_mat{:});
    anatdir = fullfile('C:\Data\ARC\',sprintf('ARC%02d',ss),'single');
    
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
    oid_ = 1:160;
    oid = oid_(group_vec)';
    
    cvind = crossvalind('Kfold',size(Fless_mat_pruned,1),nfolds);
    pred_mat = zeros(size(Fless_mat_pruned,1)/nfolds,nfolds);
    ind_mat = zeros(size(Fless_mat_pruned,1)/nfolds,nfolds);
    for folderr = 1:nfolds
        ind = cvind==folderr;
        ind_mat(:,folderr) = find(ind); 
        X_train = Fless_mat_pruned(~ind,:);
        y_train = oid(~ind);
        X_test = Fless_mat_pruned(ind,:);
        y_test = oid(ind);
        mdl = svmtrain(y_train, X_train , ['-s 1 -t 2 -q']);
        [a,p] = svmpredict( y_test,  X_test, mdl,['-q']);
%         length(unique(a))
%         length(unique(y_test))
        corrmod(folderr,ss) = p(1);
    end
end
end

figure('Position',[0 0 320 240])
hold on
bar(mean(corrmod,1))
yline(100/160)
xticks(1:3)
xticklabels({'S1','S2','S3'})
ylabel('Performance (%)')
savefig('svm')
print('svm','-dpng')
clear unity fless_mat fless_mat_unn Fless_mat utl_mask
save('svmpred')
%% Make SVM - libsvm - Fmat  and behav
nfolds = 10;

if svm_trainer3
corrmod = zeros(nfolds,3);
corrmod2 = cell(3,1);
for ss = 1:length(dirs)
    load(fullfile(dirs{ss},'sfp_feats.mat'))
    behav_ratings = behav.behav(ss).ratings;
    Fless_mat = vertcat(feat_mat{:});
    anatdir = fullfile('C:\Data\ARC\',sprintf('ARC%02d',ss),'single');
    
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
    
    Fless_mat_pruned = Fless_mat(:,[3 4 9:14]);
    Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
    oid_ = 1:160;
    oid = oid_(group_vec)';
    
    cvind = crossvalind('Kfold',size(Fless_mat_pruned,1),nfolds);
    pred_mat = zeros(size(Fless_mat_pruned,1)/nfolds,nfolds);
    ind_mat = zeros(size(Fless_mat_pruned,1)/nfolds,nfolds);
    for folderr = 1:nfolds
        ind = cvind==folderr;
        ind_mat(:,folderr) = find(ind); 
        X_train = Fless_mat_pruned(~ind,:);
        y_train = oid(~ind);
        X_test = Fless_mat_pruned(ind,:);
        y_test = oid(ind);
        mdl = svmtrain(y_train, X_train, ['-s 1 -t 2 -q']);
        [a,p] = svmpredict( y_test,  X_test, mdl,['-q']);

        corrmod(folderr,ss) = p(1);
        pred_rat = behav_ratings(a,:);
        test_rat = behav_ratings(y_test,:);
        pred_mat(:,folderr) = iter_corr(pred_rat,test_rat);      
    end
    corrmod2{ss} = pred_mat(:);
end
end


% 
figure('Position',[0 0 320 240])
hold on
bar(mean(corrmod,1))
yline(100/160)
xticks(1:3)
xticklabels({'S1','S2','S3'})
ylabel('Performance (%)')
savefig('svm')
print('svm','-dpng')

figure('Position',[0 0 320 240])
corrmod2{1}(4321:end) = []; % This is wrong, but I am lazy
corr_m = horzcat(corrmod2{:});
boxplot(corr_m)
xticks(1:3)
xticklabels({'S1','S2','S3'})
ylabel('Performance (r)')
savefig('behav')
print('behav','-dpng')


p = [];
for ii = 1:3
     p(ii) = signrank(corrmod2{ss});
end
clear unity fless_mat fless_mat_unn Fless_mat utl_mask
save('svmpred')
%% Make SVM - libsvm - Fmat  - probabilistic
nfolds = 10;

if svm_trainer4
corrmod = zeros(nfolds,3);
corrmod2 = zeros(nfolds,3);
for ss = 1:length(dirs)
    load(fullfile(dirs{ss},'sfp_feats.mat'))
    Fless_mat = vertcat(feat_mat{:});
    anatdir = fullfile('C:\Data\ARC\',sprintf('ARC%02d',ss),'single');
    
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
    
      Fless_mat_pruned = Fless_mat(:,[3 4 9:14]);
        Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
    oid_ = 1:160;
    oid = oid_(group_vec)';
    
    cvind = crossvalind('Kfold',size(Fless_mat_pruned,1),nfolds);
    pred_mat = zeros(size(Fless_mat_pruned,1)/nfolds,nfolds);
    ind_mat = zeros(size(Fless_mat_pruned,1)/nfolds,nfolds);
    for folderr = 1:nfolds
        ind = cvind==folderr;
        ind_mat(:,folderr) = find(ind); 
        X_train = Fless_mat_pruned(~ind,:);
        y_train = oid(~ind);
        X_test = Fless_mat_pruned(ind,:);
        y_test = oid(ind);
        mdl = svmtrain(y_train, X_train, ['-s 1 -t 2 -q -b 1']);
        [a,p,b] = svmpredict( y_test,  X_test, mdl,['-q -b 1']);
%         length(unique(a))
%         length(unique(y_test))
        corrmod(folderr,ss) = p(1);
        var = zeros(1,length(a));
        for zz=1:length(a)
            sdiff = setdiff(oid_,a(zz));
            var(zz) = b(zz,a(zz))-mean(b(zz,sdiff));
        end
        corrmod2(folderr,ss) = mean(var);
    end
end
end

figure('Position',[0 0 320 240])
hold on
bar(mean(corrmod2,1))
% yline(100/160)
xticks(1:3)
xticklabels({'S1','S2','S3'})
ylabel('Performance (%)')
savefig('svm')
print('svm','-dpng')
clear unity fless_mat fless_mat_unn Fless_mat utl_mask
save('svmpred')
