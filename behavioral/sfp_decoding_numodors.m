mainroot = 'C:\Work\SFP\Scripts';
rootf = fullfile(mainroot,'supporting_files');
savepath = fullfile(mainroot,'\examples\example_SVM');
nodor = 160;
wind = 7500; % Number of samples
num_subjects = 3;
behav = load(fullfile(rootf,'NEMO_perceptual2.mat'));

niter = (5:5:160)';
ntrials = length(niter);

%% SVM Decoding
nfolds = 10;
numpcs = [13 11 11]; % 90% Variance

% Additional analyses on odor distribution
corrmod = zeros(3,ntrials);
corrmodp = zeros(3,ntrials);
rng(999)

for ss = 1:num_subjects
    fprintf('subject %02d\n',ss)
    subdir = fullfile(rootf,sprintf('sfp_behav_s%02d_correct',ss));
    load(fullfile( subdir,'sfp_feats_main.mat'))

    if ss==3; s2 = 4; else; s2 = ss; end
    onsets_mat = load(fullfile( subdir,sprintf('conditions_NEMO%02d.mat',s2)),'onsets','names');
    onsets = onsets_mat.onsets;
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
    Fless_mat_pruned  = Fless_mat(:,[3 4 9:21 23:31]);
    % Fless_mat_pruned = Fless_mat(:,[3 4 9:31]);

    Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
    oid_ = 1:nodor;
    oid = oid_(group_vec)';
    Fless_mat_pruned = zscore(Fless_mat_pruned,1);

    for tt = 1:ntrials
        oid_select = sort(randperm(160,niter(tt)));
        trial_cutoff = ismember(oid,oid_select);

        [coeff,Fless_mat_pruned_thisrun,~,~,var] = pca(Fless_mat_pruned(trial_cutoff,:));
        ncomp = sum(cumsum(var)<90)+1;
        Fless_mat_pruned_thisrun = Fless_mat_pruned_thisrun(:,1:ncomp);

        [~,predictions_vec] = Classify_Permute_VS2(Fless_mat_pruned_thisrun, oid( trial_cutoff ), 5);
        % behavioral_corr
        accuracies = predictions_vec==oid( trial_cutoff );

        % Accuracies
        corrmod(ss,tt) = sum(accuracies)/length(accuracies)-1/niter(tt);
        corrmodp(ss,tt) = ARC_computePValueOneTailed(corrmod(ss,tt), niter(tt), length(accuracies));
    end
end

figure()
hold on
plot(niter,corrmod')

%% Odor selection - training
ss = 2;
subdir = fullfile(rootf,sprintf('sfp_behav_s%02d_correct',ss));
load(fullfile( subdir,'sfp_feats_main.mat'))
Fless_mat = vertcat(feat_mat{:});
Fless_mat_pruned  = Fless_mat(:,[3 4 9:21 23:31]);
Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
Fless_mat_pruned = zscore(Fless_mat_pruned,1);
% 
% [coeff,Fless_mat_pruned,~,~,var] = pca(Fless_mat_pruned);
% Fless_mat_pruned = Fless_mat_pruned(:,1:numpcs(ss));
    
% Additional analyses on odor distribution
niter = 10*ones(1,1000);
ntrials = length(niter);
corrmod = zeros(3,ntrials);
corrmodp = zeros(3,ntrials);
corr_oids = zeros(10,ntrials);
rng(999)
for tt = 1:ntrials
    if mod(tt,20)==0; fprintf('Running iter: %03d\n',tt); end
    oid_select = sort(randperm(160,niter(tt)));
    trial_cutoff = ismember(oid,oid_select);
    corr_oids(:,tt) = oid_select;


    [coeff,Fless_mat_pruned_thisrun,~,~,var] = pca(Fless_mat_pruned(trial_cutoff,:));
    ncomp = sum(cumsum(var)<90)+1;
    Fless_mat_pruned_thisrun = Fless_mat_pruned_thisrun(:,1:ncomp);

    [~,predictions_vec] = Classify_Permute_VS2(Fless_mat_pruned_thisrun, oid( trial_cutoff ), 5);
    % behavioral_corr
    accuracies = predictions_vec==oid( trial_cutoff );

    % Accuracies
    corrmod(ss,tt) = sum(accuracies)/length(accuracies)-1/niter(tt);
    corrmodp(ss,tt) = ARC_computePValueOneTailed(corrmod(ss,tt), niter(tt), length(accuracies));
end

[~, sort_id] = sort( corrmod(ss,:),'descend');
corr_oids_sort = corr_oids(:,sort_id);

%% Odor selection - validation

ss = 2;
subdir = fullfile(rootf,sprintf('sfp_behav_s%02d_correct',ss));
load(fullfile( subdir,'sfp_feats_main.mat'))
Fless_mat = vertcat(feat_mat{:});
Fless_mat_pruned  = Fless_mat(:,[3 4 9:21 23:31]);
Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
Fless_mat_pruned = zscore(Fless_mat_pruned,1);
% 
% [coeff,Fless_mat_pruned,~,~,var] = pca(Fless_mat_pruned);
% Fless_mat_pruned = Fless_mat_pruned(:,1:numpcs(ss));
    
% Additional analyses on odor distribution

rng(999)
% 
% oid_select_l = [240 8785 7335 7410 15510 7720 1068 9862 10364 5365027];
oid_select_l = [240 8785 7335 7410 10430 7720 1001 9862 10364 5365027];
M = containers.Map(behav(2).cid,1:160);
oid_select = []; for zz = 1:10; oid_select(zz)=M(oid_select_l(zz)); end


trial_cutoff = ismember(oid,oid_select);

[coeff,Fless_mat_pruned_thisrun,~,~,var] = pca(Fless_mat_pruned(trial_cutoff,:));
ncomp = sum(cumsum(var)<90)+1;
Fless_mat_pruned_thisrun = Fless_mat_pruned_thisrun(:,1:ncomp);

[~,predictions_vec] = Classify_Permute_VS2(Fless_mat_pruned_thisrun, oid( trial_cutoff ), 5);
% behavioral_corr
accuracies = predictions_vec==oid( trial_cutoff );

% Accuracies
corrmod = sum(accuracies)/length(accuracies);
corrmodp = ARC_computePValueOneTailed(corrmod, 10, length(accuracies));




