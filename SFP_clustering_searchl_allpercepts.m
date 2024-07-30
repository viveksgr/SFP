%% General Settings
root = 'C:\Work';
tic
settings_.nodor = 160;
settings_.wind = 7500; % Number of samples

% m_id = [8 13; 8 15; 11 16]; % N-settings_.wind
% m_id = [13 61; 15 42; 16 41]; % N-settings_.wind
% m_id = [3 4; 3 4; 3 4]; % Inhale Exhale
settings_.ncomps = [5 5 5 5];
settings_.nsniffcomp = 31;
settings_.featorflessnot = false;
settings_.numClusters = 60;
settings_.pcamaker = false;
settings_.mapper = true;

dirs = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
    fullfile(root ,'\SFP\sfp_behav_s02_correct');
    fullfile(root ,'\SFP\sfp_behav_s04_correct')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
    fullfile(root,'ARC\ARC\ARC02\single');
    fullfile(root,'ARC\ARC\ARC03\single')};

savepath = 'C:\Work\SFP\Clustering\allpercepts';
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3.nii';

anat_names = {'PC','AMY','OFC','OT','AON'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwOT.nii','rwAON.nii'};
nanat = length(anat_names);

settings_.fmasker  = true; % Turn on the functional mask
settings_.single_n = false; % Noisepool
settings_.single_c = true; % Cutoff from sign voxels
settings_.mapper   = true;
settings_.loadvec  = [3 4 9:21 23:settings_.nsniffcomp];

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')
rsa_P1_ = cell(3,nanat,3); % Complete
hold on
anat_cell = {};
% Subject - index
behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));

hold on
sbplt = 0;
numpcs = [13 11 11]; % 90% Variance
idxmats = cell(3,2);
nfeatures = length(settings_.loadvec);
ndescrip = 18;
wt_sc_mat = zeros(nfeatures,ndescrip,nanat,3);
tt_sc_mat = zeros(nfeatures,ndescrip,nanat,3);

for ss = [1 2 3] % Subject
    fprintf('Subject: %02d\n',ss)
    h = waitbar(0, 'Please wait...'); % Progress bar per subject
    if ss==3; s2 = 4; else; s2 = ss; end
    statpath = dirs{ss};
    anatdir = dirs2{ss};
    % savepath = dirs3{ss};
    mkdir(savepath)
  
    load(fullfile(statpath,'sfp_feats_main.mat'))
    load(fullfile(statpath,'task_struct_trialwise.mat'))
    Fless_mat = vertcat(fless_mat{:});
    Fless_mat_pruned = Fless_mat(:,1:settings_.wind);

    Feat_mat_pruned = vertcat(feat_mat{:});
    Feat_mat_pruned =  Feat_mat_pruned(:,[settings_.loadvec]) ;

    if settings_.pcamaker
        Feat_mat_pruned(isnan(Feat_mat_pruned))=0;
       Feat_mat_pruned = zscore(Feat_mat_pruned,1);

        [coeff,sc,~,~,var] = pca(Feat_mat_pruned);
        Feat_mat_pruned = sc(:,1:numpcs(ss));
    end

    %% Gray Matter, Functional and Anatomical Masks
    mask = (spm_read_vols(spm_vol(fullfile(anatdir, maskfile)))); % Mask used to construct odor files
    mask(isnan(mask))=0;
    mask = logical(mask);
    fmask = (spm_read_vols(spm_vol(fullfile(anatdir, fmaskfile)))); % Mask used to examine voxels in RSA
    fmask(isnan(fmask))=0;
    if ~settings_.fmasker
        fmask = fmask +0.1;
    end
    fmask = logical(fmask); % Only choose voxels with significant odor evoked activity
    fmask_1d = fmask(mask);
    marea = and(fmask,mask);
    anatpath = anatdir;

    % Model names
    masks_set = [];
    masks_set_cell = {};

    anatmasks = [];
    for ii = 1:length(anat_masks)
        m1 = spm_read_vols(spm_vol(fullfile(anatpath,anat_masks{ii})));
        m1(isnan(m1))=0;
        m1(m1<=0.01)=0;
        m1(m1>0) = 1;
        m1 = m1(mask);
        masks_set(:,ii)=m1(fmask_1d);

        masks_set_cell{ii} = fmask_1d(logical(m1));
        anatmask = (spm_read_vols(spm_vol(fullfile(anatpath, anat_masks{ii}))));
        anatmask(isnan(anatmask)) = 0;
        anatmask = logical(anatmask);
        anatmask = and(anatmask,marea);
        anatmasks(:,:,:,ii) = anatmask;
    end
    anat_cell{ss}= anatmasks;
    masks_set(isnan(masks_set))=0;
    tnvox = sum(anatmasks,'all');

    %% Defining RDMs
    if settings_.featorflessnot; mainmat = Feat_mat_pruned; dister = 'sqeuclidean'; else; mainmat = Fless_mat_pruned(:,1:4000); dister = 'correlation';  end
    mainmat(isnan(mainmat))=0;
    if settings_.featorflessnot
        mainmat = zscore(mainmat,[],1);
    end
  
    fname_ = 'conditions_NEMO%02d.mat'; 
    onsets = load(fullfile(anatdir,sprintf(fname_,s2)),'onsets');
    onsets = onsets.onsets;

    onsets = load(fullfile(anatdir,sprintf('conditions_NEMO%02d.mat',s2)),'onsets');
    onsets = onsets.onsets;
    group_vec = cell(settings_.nodor,1);
    unity = [];
    for ii2 = 1:settings_.nodor
        group_vec{ii2} = ii2*ones(length(onsets{ii2}),1);
        unity = blkdiag(unity,ones(length(onsets{ii2})));
    end
    group_vec = vertcat(group_vec{:});
    [onsetlist,argsort] = sort(vertcat(onsets{:}));
    group_vec = group_vec(argsort);
    unity = unity(argsort,argsort);

    % Behavioral RSMs
    behav_ratings = behav.behav(ss).ratings;
    behav_ratings = behav_ratings(group_vec,:);

    ndescrip = size(behav_ratings,2);


    for pp = 1:ndescrip

        [idx_int,~] = kmeans(behav_ratings(:,pp),settings_.numClusters,'MaxIter',1000,'Replicates',1000);
        Aint = SFP_splitapply_mean(mainmat,idx_int);
        utl_mask = logical(triu(ones(length(unique(idx_int))),1)); % All possible odor
        Aint_corr = SFP_extractUpperTriangles(   Aint,utl_mask);

        % Nuisance regressors
        task_run1 = SFP_splitapply_mean(SFP_splitapply_mean(task_run',idx_int)',idx_int);
        sess_run1 = SFP_splitapply_mean(SFP_splitapply_mean(sess_run',idx_int)',idx_int);
        set_run1 = SFP_splitapply_mean(SFP_splitapply_mean(set_run',idx_int)',idx_int);
     
        %% Representational connectivity
        kvox = 0;
        map_area = zeros([size(anat_cell{ss}) 3]);
        for ii = 1:length(anat_names)
        fprintf('area:%02d\n',ii)
        modelmd_ = load(fullfile(anatdir,'desniff',anat_names{ii},'TYPEC_FITHRF_GLMDENOISE.mat'),'modelmd','noisepool');
        modelmd = squeeze(modelmd_.modelmd);
        noisepool = modelmd_.noisepool;
        if settings_.single_c
            modelmd = modelmd(masks_set_cell{ii},:);
            noisepool = noisepool(masks_set_cell{ii});
        end
        S_omat_vals_r = modelmd;
        [r1,~] = find(isnan(S_omat_vals_r));
        S_omat_vals_r(r1,:) = [];
        % M_anat = 1-corrcoef(S1_omat_vals);

        M1 = SFP_splitapply_mean(S_omat_vals_r',idx_int);
        M1 = zscore(M1,[],1);
        M1_anat = corrcoef(M1');

        [wt,t_sc1] = ARC_multicomputeWeights_tsc([Aint_corr task_run1(utl_mask) sess_run1(utl_mask) ], M1_anat(utl_mask));
        wt_sc_mat(:,pp,ii,ss)=wt(2:end-2);
        tt_sc_mat(:,pp,ii,ss)=t_sc1(2:end-2);
        end
        
    end
end

rsa_P1 = squeeze(sum(abs(tt_sc_mat)>tinv(0.975,nchoosek(60,2))));
rsa_P1_mat = cat(1,rsa_P1(1:2,:,:),mean(rsa_P1(3:end,:,:)));

ARC_barplot(permute(rsa_P1_mat,[3 2 1]),true)
xticks(1:5)
xticklabels(anat_names)
legend({'Int','Pls','Percepts'})
ylabel('Num sniff descripts')
savefig('numfeat')
print('numfeat','-dpng')

% Comparision with behavioral tasks
load('C:\Work\SFP\Final_plots\Behavioral\Multilogistic\nominal\dataset.mat','wt_mat_cell')

wt_corr = zeros(3,nanat,3);
for ss = 1:3
    for ii = 1:nanat
        wt_corr(ss,ii,1)= fastcorr(abs(squeeze(wt_sc_mat(:,1,ii,ss))),abs(wt_mat_cell{ss}(:,1)));
        wt_corr(ss,ii,2)= fastcorr(abs(squeeze(wt_sc_mat(:,2,ii,ss))),abs(wt_mat_cell{ss}(:,2)));
        wt_corr(ss,ii,3)= fastcorr(abs(squeeze(mean(wt_sc_mat(:,3:end,ii,ss)))),abs(mean(wt_mat_cell{ss}(:,3:end))));
    end
end
ARC_barplot(wt_corr)
xticks(1:5)
xticklabels(anat_names)
legend({'Int','Pls','Percepts'})
ylabel('Similarity with behavioral model')
savefig('corrfeatabs')
print('corrfeatabs','-dpng')
