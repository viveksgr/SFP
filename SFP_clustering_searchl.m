%% General Settings
root = 'C:\Work';
tic
settings_.nodor = 160;
settings_.wind = 7500; % Number of samples

% m_id = [8 13; 8 15; 11 16]; % N-settings_.wind
% m_id = [13 61; 15 42; 16 41]; % N-settings_.wind
% m_id = [3 4; 3 4; 3 4]; % Inhale Exhale
settings_.ncomps = [5 5 5 5];
settings_.nsniffcomp = 32;
settings_.featorflessnot = true;
settings_.numClusters = 60;
settings_.numClusters2 = 60;
settings_.clustbehav = true;
settings_.sniffupdate = false;
settings_.mapper = false;

dirs = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
    fullfile(root ,'\SFP\sfp_behav_s02_correct');
    fullfile(root ,'\SFP\sfp_behav_s04_correct')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
    fullfile(root,'ARC\ARC\ARC02\single');
    fullfile(root,'ARC\ARC\ARC03\single')};

savepath = 'C:\Work\SFP\Clustering\Feat_main_updated2';
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
rsa_P1_ = cell(3,nanat,2); % Complete
hold on
anat_cell = {};
% Subject - index
behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));

hold on
sbplt = 0;
for ss = [1 2 3] % Subject
    fprintf('Subject: %02d\n',ss)
    h = waitbar(0, 'Please wait...'); % Progress bar per subject
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
        [idx,A1] = kmeans(mainmat,settings_.numClusters,'Distance',dister,'Replicates',1);
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

    [idx2,A1] = kmeans(behav_ratings,settings_.numClusters2,'Replicates',1);
    if settings_.clustbehav; group_vec = idx2; end
    A2 = SFP_splitapply_mean(mainmat,group_vec);
    A2_corr = corrcoef(A2');

    % Nuisance regressors
    task_run1 = SFP_splitapply_mean(SFP_splitapply_mean(task_run',idx)',idx);
    sess_run1 = SFP_splitapply_mean(SFP_splitapply_mean(sess_run',idx)',idx);
    set_run1 = SFP_splitapply_mean(SFP_splitapply_mean(set_run',idx)',idx);
    task_run2 = SFP_splitapply_mean(SFP_splitapply_mean(task_run',group_vec)',group_vec);
    sess_run2 = SFP_splitapply_mean(SFP_splitapply_mean(sess_run',group_vec)',group_vec);
    set_run2 = SFP_splitapply_mean(SFP_splitapply_mean(set_run',group_vec)',group_vec);


    %% Representational connectivity
    kvox = 0;
    for ii = 1:length(anat_names)
        fprintf('area:%02d\n',ii)
        modelmd_ = load(fullfile(anatdir,'desniff',anat_names{ii},'TYPEC_FITHRF_GLMDENOISE.mat'),'modelmd','noisepool');
        modelmd = squeeze(modelmd_.modelmd);
        noisepool = modelmd_.noisepool;
        if settings_.single_c
            modelmd = modelmd(masks_set_cell{ii},:);
            noisepool = noisepool(masks_set_cell{ii});
        end
        S_omat_vals = modelmd;
        [r1,~] = find(isnan(S_omat_vals));
        S_omat_vals(r1,:) = [];
        % M_anat = 1-corrcoef(S1_omat_vals);

        rget = 3;
        sz = size(logical(anatmasks(:,:,:,ii)));
        ref_vox = round(sz/2);
        [MX, MY, MZ] = ndgrid(1:sz(1),1:sz(2),1:sz(3));
        radii = sqrt((MX-ref_vox(1)).^2 + (MY-ref_vox(2)).^2 + (MZ-ref_vox(3)).^2);
        % prototype sphere index for radii<rget that can be used everywhere in the brain
        radius_index = find(radii<rget) - sub2ind(sz,ref_vox(1),ref_vox(2),ref_vox(3)); % Neighbors index

        % Gray mask configuration
        lin_index = find(logical(anatmasks(:,:,:,ii)));
        nvox=length(lin_index);
        % resultsvolume index
        linindexconv = zeros(sz);
        linindexconv(lin_index) = 1:length(lin_index);

        rsa_vec_1 = zeros(nvox,1);
        rsa_vec_2 = zeros(nvox,1);

        for cnt2 = 1:nvox

            indexindex2 = radius_index + lin_index(cnt2);
            indexindex2 = intersect(lin_index, indexindex2);
            indexindex2 = linindexconv(indexindex2);
            S_omat_vals_r = S_omat_vals(indexindex2,:);
            [r1,~] = find(isnan(S_omat_vals_r));
            S_omat_vals_r(r1,:) = [];

            if (size(S_omat_vals_r,1)>1)


                M1 = SFP_splitapply_mean(S_omat_vals_r',idx);
                M1 = zscore(M1,[],1);
                M2 = SFP_splitapply_mean(S_omat_vals_r',group_vec);
                M2 = zscore(M2,[],1);

                M1_anat = corrcoef(M1');
                M2_anat = corrcoef(M2');


                utl_mask = logical(triu(ones(length(unique(idx))),1)); % All possible odors
                [~,t_sc] = ARC_multicomputeWeights_tsc([A1_corr(utl_mask) task_run1(utl_mask) sess_run1(utl_mask) set_run1(utl_mask)], M1_anat(utl_mask));
                % [~,t_sc] = ARC_multicomputeWeights_tsc([A1_corr(utl_mask)], M1_anat(utl_mask));

                rsa_vec_1(cnt2) = t_sc(2);

                utl_mask = logical(triu(ones(length(unique(group_vec))),1)); % All possible odors
                [~,t_sc] = ARC_multicomputeWeights_tsc([A2_corr(utl_mask) task_run2(utl_mask) sess_run2(utl_mask) set_run2(utl_mask)], M2_anat(utl_mask));
                % [~,t_sc] = ARC_multicomputeWeights_tsc([A2_corr(utl_mask)], M2_anat(utl_mask));

                rsa_vec_2(cnt2) = t_sc(2);
            end
            kvox = kvox+1;
            waitbar(kvox/(tnvox), h, sprintf('Processing S:%02d, %.2f %%', ss, (kvox/(tnvox)*100))) % Update progress
        end

        rsa_P1_{ss,ii,1} = rsa_vec_1;
        rsa_P1_{ss,ii,2} = rsa_vec_2;

        if settings_.mapper
            rsa_vec_1 = unmasker(rsa_vec_1,logical(anatmasks(:,:,:,ii)));
            rsa_vec_2 = unmasker(rsa_vec_2,logical(anatmasks(:,:,:,ii)));
            write_reshaped_nifty(rsa_vec_1, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_peak_%s',ss,anat_names{ii}));
            write_reshaped_nifty(rsa_vec_2, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_RSA_%s',ss,anat_names{ii}));
            write_reshaped_nifty(rsa_vec_1-rsa_vec_2, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_peak-RSA_%s',ss,anat_names{ii}));
            write_reshaped_nifty(rsa_vec_2-rsa_vec_1, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_RSA-peak_%s',ss,anat_names{ii}));
        end

    end
    close(h)
end


rsa_P1 = cellfun(@(x) (sum(x> tinv(0.99,sum(utl_mask,'all')))./length(x))*100,rsa_P1_);
rsa_P1(3,5,:) = nan;
ARC_barplot(rsa_P1)
gcf
% yline(r2t(0.05,length(M1_anat_vec)))
xticks(1:nanat)
xticklabels(anat_names);
% r2t(0.05,length(A2_corr_vec))
legend('S_I','S_O')
ylabel('Representational Similarity (%)')
savefig(fullfile(savepath,'feat_map'))
print(fullfile(savepath,'feat_map'),'-dpng')

rsa_P1_comp = SFP_computePercentages(rsa_P1_,tinv(0.95,sum(utl_mask,'all')));
ARC_barplot(rsa_P1_comp)
gcf
% yline(r2t(0.05,length(M1_anat_vec)))
xticks(1:nanat)
xticklabels(anat_names);
% r2t(0.05,length(A2_corr_vec))
legend('S_I>S_O','S_O>S_I','Unclassified')
ylabel('Representational Similarity (%)')
savefig(fullfile(savepath,'feat_map_comp'))
print(fullfile(savepath,'feat_map_comp'),'-dpng')

clearvars fless_mat mainmat Fless_mat_pruned Fless_mat unity task_run set_run sess_run anat_cell

% clear Fmat_1_m behav_mat unity M_anat M_reg n_reg M_unstack fless_mat fless_mat_unn modelmd_ modelmd S_omat_vals utl_mask utl_mask2
save(fullfile(savepath,'ARC_RSA'),'settings_','rsa_P1')
toc
