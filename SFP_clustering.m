%% General Settings
root = 'C:\Work';
tic
settings_.nodor = 160;
settings_.wind = 7500; % Number of samples

% m_id = [8 13; 8 15; 11 16]; % N-settings_.wind
% m_id = [13 61; 15 42; 16 41]; % N-settings_.wind
% m_id = [3 4; 3 4; 3 4]; % INhale Exhale
settings_.ncomps = [5 5 5 5];
settings_.nsniffcomp = 32;
settings_.featorflessnot = true; 
settings_.numClusters = 60;
settings_.numClusters2 = 60;
settings_.clustbehav = true;
settings_.sniffupdate = false;
settings_.loadcluster = true;

dirs = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
        fullfile(root ,'\SFP\sfp_behav_s02_correct');
        fullfile(root ,'\SFP\sfp_behav_s04_correct')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
        fullfile(root,'ARC\ARC\ARC02\single');
        fullfile(root,'ARC\ARC\ARC03\single')};

savepath = 'C:\Work\SFP\Clustering\Feat_main_nonsearchl';
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3.nii'; 

anat_names = {'PC','AMY','OFC','OT','AON'}; 
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwOT.nii','rwAON.nii'};
nanat = length(anat_names);

settings_.fmasker = true; % Turn on the functional mask
settings_.single_n = false; % Noisepool
settings_.single_c = true; % Cutoff from sign voxels
settings_.mapper = true;
settings_.loadvec = [3 4 9:settings_.nsniffcomp ];

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')
rsa_P1 = zeros(3,nanat,2); % Complete
hold on
anat_cell = {};
% Subject - index
behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));

hold on
sbplt = 0;
totalMI = zeros(3,nanat);
pvalues = zeros(3,nanat);

for ss = [1 2 3] % Subject
    fprintf('Subject: %02d\n',ss)
    if ss==3; s2 = 4; else; s2 = ss; end
    statpath = dirs{ss};
    anatdir = dirs2{ss};
    % savepath = dirs3{ss};
    % mkdir(savepath)

    load(fullfile(statpath,'sfp_feats_main.mat'))
    load(fullfile(statpath,'task_struct_trialwise.mat'))
    Fless_mat = vertcat(fless_mat{:});
    Fless_mat_pruned = Fless_mat(:,1:settings_.wind);
    
    Feat_mat_pruned = vertcat(feat_mat{:});
    Feat_mat_pruned =  Feat_mat_pruned(:,[settings_.loadvec]) ;

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
     if settings_.featorflessnot; mainmat = Feat_mat_pruned; dister = 'sqeuclidean'; else; mainmat = Fless_mat_pruned; dister = 'sqeuclidean';  end
    mainmat(isnan(mainmat))=0;
    if settings_.featorflessnot
        mainmat = zscore(mainmat,[],1);
    end

    if settings_.sniffupdate; fname_ = 'conditions_sniffupdate_NEMO%02d.mat'; else; fname_ = 'conditions_NEMO%02d.mat'; end
    onsets = load(fullfile(anatdir,sprintf(fname_,s2)),'onsets');
    onsets = onsets.onsets;

    if settings_.sniffupdate
        idx =  SFP_mapMainToSniff(onsets);
        A1 = SFP_splitapply_mean(mainmat,idx);
    else
        if settings_.loadcluster
            load(fullfile('C:\Work\SFP\Clustering\Feat_main_updated_tsc','ARC_RSA.mat'),'idxmats')
            idx = idxmats{ss,1};
            A1 =  SFP_splitapply_mean(mainmat,idx);
        else 
            [idx,A1] = kmeans(mainmat,settings_.numClusters,'Distance',dister,'Replicates',1);
        end
    end
    A1_corr = corrcoef(A1');

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
    if settings_.loadcluster
        idx2 = idxmats{ss,2};
    else
        [idx2] = kmeans(behav_ratings,settings_.numClusters2,'Replicates',1);
    end
    if settings_.clustbehav; group_vec = idx2; end
    A2 = SFP_splitapply_mean(mainmat,group_vec);
    A2_corr = corrcoef(A2');

    %% Representational connectivity
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

        M1 = SFP_splitapply_mean(S_omat_vals_r',idx);
        M1 = zscore(M1,[],1);
        M2 = SFP_splitapply_mean(S_omat_vals_r',group_vec);
        M2 = zscore(M2,[],1);
        [totalMI(ss,ii), pValues(ss,ii)] = sfP_MI(M1, M2, 1000);



        M1_anat = corrcoef(M1');
        M2_anat = corrcoef(M2');
       
        % Nuisance regressors
        task_run1 = SFP_splitapply_mean(SFP_splitapply_mean(task_run',idx)',idx);
        sess_run1 = SFP_splitapply_mean(SFP_splitapply_mean(sess_run',idx)',idx);
        set_run1 = SFP_splitapply_mean(SFP_splitapply_mean(set_run',idx)',idx);
        task_run2 = SFP_splitapply_mean(SFP_splitapply_mean(task_run',group_vec)',group_vec);
        sess_run2 = SFP_splitapply_mean(SFP_splitapply_mean(sess_run',group_vec)',group_vec);
        set_run2 = SFP_splitapply_mean(SFP_splitapply_mean(set_run',group_vec)',group_vec);   
        
        utl_mask = logical(triu(ones(length(unique(idx))),1)); % All possible odors       
        [~,t_sc] = ARC_multicomputeWeights_tsc([A1_corr(utl_mask) task_run1(utl_mask) sess_run1(utl_mask) set_run1(utl_mask)], M1_anat(utl_mask));
        % [~,t_sc] = ARC_multicomputeWeights_tsc([A1_corr(utl_mask)], M1_anat(utl_mask));
       
        rsa_P1(ss,ii,1) = t_sc(2);
             
        utl_mask = logical(triu(ones(length(unique(group_vec))),1)); % All possible odors       
        [~,t_sc] = ARC_multicomputeWeights_tsc([A2_corr(utl_mask) task_run2(utl_mask) sess_run2(utl_mask) set_run2(utl_mask)], M2_anat(utl_mask));
         % [~,t_sc] = ARC_multicomputeWeights_tsc([A2_corr(utl_mask)], M2_anat(utl_mask));
       
        rsa_P1(ss,ii,2) = t_sc(2);
    
    end
end

rsa_pvals = tcdf(rsa_P1,sum(utl_mask,'all'),'upper');

S_mat = squeeze(mean(rsa_P1));
rsa_pvals_mean = tcdf(S_mat,sum(utl_mask,'all'),'upper');

S_err = squeeze(std(rsa_P1))./sqrt(3);
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
yline(tinv(0.95,sum(utl_mask,'all')))
xticks(1:nanat)
xticklabels(anat_names);
% legend({'Perceptual','Chemical','Mutual'})
% legend()
% Subject data points
c_s = {'r','g','b'}; % Data dots for subjects
for ii = 1:nanat % For bars for perceptual, chemical and combinations
    for jj = 1:3
        plot(x_m(:,ii),squeeze(rsa_P1(jj,ii,:)),c_s{jj},'handle','off')
    end
end
legend('S_I','S_O')
ylabel('Representational Similarity (t)')
% yline(r2t(0.05,sum(utl_mask2(:))));
% yline(r2t(0.05,nchoosek(length( group_vec),2)));
mkdir(savepath)

savefig(fullfile(savepath,'feat_map'))
print(fullfile(savepath,'feat_map'),'-dpng')
clearvars fless_mat mainmat Fless_mat_pruned Fless_mat unity task_run set_run sess_run anat_cell

% clear Fmat_1_m behav_mat unity M_anat M_reg n_reg M_unstack fless_mat fless_mat_unn modelmd_ modelmd S_omat_vals utl_mask utl_mask2
save(fullfile(savepath,'ARC_RSA'),'settings_','rsa_P1') 
toc
%% Kmeans elbow method
elbow =false;
if elbow
figure()
hold on
krange = 5:5:160;
for ss = [1 2 3] % Subject
    disp('Running k-means elbow')
    fprintf('Subject: %02d\n',ss)
    if ss==3; s2 = 4; else; s2 = ss; end
    statpath = dirs{ss};
    anatdir = dirs2{ss};
    % savepath = dirs3{ss};
    % mkdir(savepath)

    load(fullfile(statpath,'sfp_feats_main.mat'))
    Fless_mat = vertcat(fless_mat{:});
    Fless_mat_pruned = Fless_mat(:,1:settings_.wind);
   
    wcss = SFP_kmeansElbowMethod(Fless_mat_pruned, krange ,'correlation');
    subplot(3,3,ss)
    plot(krange,wcss)
    title(sprintf('Sub:%02d, Fless',ss))

    Feat_mat_pruned = vertcat(feat_mat{:});
    Feat_mat_pruned =  Feat_mat_pruned(:,[settings_.loadvec]) ;
    wcss = SFP_kmeansElbowMethod(Feat_mat_pruned, krange );
    
    subplot(3,3,ss+3)
    plot(krange,wcss)
    title(sprintf('Sub:%02d, Fmat',ss))

    behav_ratings = behav.behav(ss).ratings;
    % behav_ratings = behav_ratings(group_vec,:);

    wcss = SFP_kmeansElbowMethod(behav_ratings, krange );
    
    subplot(3,3,ss+6)
    plot(krange,wcss)
    title(sprintf('Sub:%02d, Behavioral',ss))

end
end