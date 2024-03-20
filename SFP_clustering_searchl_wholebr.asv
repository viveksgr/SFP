%% General Settings
root = 'C:\Work';
settings_.nodor = 160;
settings_.wind = 7500; % Number of samples

% m_id = [8 13; 8 15; 11 16]; % N-settings_.wind
% m_id = [13 61; 15 42; 16 41]; % N-settings_.wind
% m_id = [3 4; 3 4; 3 4]; % INhale Exhale
settings_.ncomps = [5 5 5 5];
settings_.nsniffcomp = 14;
settings_.ncomp_plt = 8;
settings_.prfm_pca = false;
settings_.reg_mat = true;
settings_.featorflessnot = false;
settings_.numClusters = 160;
settings_.loadkmeansidx = true;

dirs = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
    fullfile(root ,'\SFP\sfp_behav_s02_correct');
    fullfile(root ,'\SFP\sfp_behav_s04_correct')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
    fullfile(root,'ARC\ARC\ARC02\single');
    fullfile(root,'ARC\ARC\ARC03\single')};

dirs3 = {fullfile(root ,'SFP\sfp_behav_s01_clusterwb');
    fullfile(root ,'\SFP\results\Clustering\sfp_behav_s02_cluster');
    fullfile(root ,'\SFP\results\Clustering\sfp_behav_s04_cluster')};

maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3.nii';

anat_names = {'PC','AMY','OFC','OT','AON'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwOT.nii','rwAON.nii'};
nanat = length(anat_names);

settings_.fmasker = true; % Turn on the functional mask
settings_.single_n = false; % Noisepool
settings_.single_c = true; % Cutoff from sign voxels
settings_.mapper = true;
settings_.loadvec = [3 4 9:14];

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')
rsa_P1_ = cell(3,nanat,2); % Complete
hold on
anat_cell = {};
% Subject - index
behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));

hold on
sbplt = 0;

for ss = [1] % Subject
    h = waitbar(0, 'Please wait...'); % Progress bar per subject
    fprintf('Subject: %02d\n',ss)
    if ss==3; s2 = 4; else; s2 = ss; end
    statpath = dirs{ss};
    anatdir = dirs2{ss};
    savepath = dirs3{ss};
    mkdir(savepath)

    load(fullfile(statpath,'sfp_feats_corrected.mat'))
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
    utl_mask = logical(triu(ones(length(unique(group_vec))),1)); % All possible odors

    onsets_sniff = load(fullfile(dirs3{ss},sprintf('conditions_sniff_NEMO%02d.mat',s2)),'onsets');
    onsets_sniff = onsets_sniff.onsets;
    % Behavioral RSMs
    behav_ratings = behav.behav(ss).ratings;
    behav_corr = corrcoef(behav_ratings');

    if settings_.featorflessnot; mainmat = Feat_mat_pruned; else; mainmat = Fless_mat_pruned; end
    mainmat(isnan(mainmat))=0;


    if settings_.featorflessnot
        mainmat = zscore(mainmat,[],1);
    end

    if settings_.loadkmeansidx
         [idx] = SFP_mapMainToSniff(onsetlist, onsets_sniff);
         A1 = SFP_splitapply_mean(mainmat,idx);
    else
        [idx,A1] = kmeans(mainmat,settings_.numClusters);
    end
    A2 = splitapply(@mean,mainmat,group_vec);
    A1_corr = corrcoef(A1');
    A2_corr = corrcoef(A2');

    %% Representational connectivity
    kvox = 0;
   
    sniffdata = load(fullfile(dirs3{ss},'full_zscored_sniffmat.mat'),'odor_responses_nn');
    sniffdata = squeeze(sniffdata.odor_responses_nn);

    odordata = load(fullfile(dirs3{ss},'full_SFP.mat'),'odor_responses_nn');
    odordata =  squeeze(odordata.odor_responses_nn(:,6,:));

    assert(size(odordata,1)==size(sniffdata,1))


    S_omat_vals_sniff = zscore(sniffdata,[],2);
    [r1,~] = find(isnan(S_omat_vals_sniff));
    S_omat_vals_odor = zscore(odordata,[],2);
    [r2,~] = find(isnan(S_omat_vals_odor));


    S_omat_vals_sniff(unique([r1 r2]),:) = [];
    S_omat_vals_odor(unique([r1 r2]),:) = [];

    nanmask = false(size(odordata,1),1);
    nanmask(unique([r1 r2])) = true;
    nanmask = logical(unmasker(nanmask,mask));
    mainmask = and(~nanmask,mask);

    assert(sum(mainmask(:))==length(S_omat_vals_odor))

    % Searchlight initiation
    rget = 3;
    sz = size(logical(mainmask));
    ref_vox = round(sz/2);
    [MX, MY, MZ] = ndgrid(1:sz(1),1:sz(2),1:sz(3));
    radii = sqrt((MX-ref_vox(1)).^2 + (MY-ref_vox(2)).^2 + (MZ-ref_vox(3)).^2);
    % prototype sphere index for radii<rget that can be used everywhere in the brain
    radius_index = find(radii<rget) - sub2ind(sz,ref_vox(1),ref_vox(2),ref_vox(3)); % Neighbors index

    % Gray mask configuration
    lin_index = find(logical(mainmask));
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

        S_omat_vals_sr = S_omat_vals_sniff(indexindex2,:);
        S_omat_vals_or = S_omat_vals_odor(indexindex2,:);


        if (size(S_omat_vals_sr,1)>1)

            % S_I estimation
            M1_anat = corrcoef( S_omat_vals_sr);
            M1_anat_vec = M1_anat(utl_mask);
            A1_corr_vec = A1_corr(utl_mask);
            [r1,~] = find(isnan(M1_anat_vec));
            M1_anat_vec(r1,:) = [];
            A1_corr_vec(r1,:) = [];
            rsa_vec_1(cnt2) = fastcorr( M1_anat_vec,A1_corr_vec);
            fastcorr(behav_corr(utl_mask),A1_corr(utl_mask))

            % S_0 estimation
            M2_anat = corrcoef( S_omat_vals_or);
            M2_anat_vec = M2_anat(utl_mask);
            A2_corr_vec = A2_corr(utl_mask);
            [r1,~] = find(isnan(M2_anat_vec));
            M2_anat_vec(r1,:) = [];
            A2_corr_vec(r1,:) = [];
            rsa_vec_2(cnt2) = fastcorr( M2_anat_vec,A2_corr_vec);
            fastcorr(behav_corr(utl_mask),A2_corr(utl_mask))

        end
        kvox = kvox+1;
        waitbar(kvox/(nvox), h, sprintf('Processing S:%02d, %.2f %%', ss, (kvox/(nvox)*100))) % Update progress
    end

    if settings_.mapper
        rsa_vec_1 = unmasker(rsa_vec_1,logical(mainmask));
        rsa_vec_2 = unmasker(rsa_vec_2,logical(mainmask));
        write_reshaped_nifty(rsa_vec_1, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_peak_%s',ss,anat_names{ii}));
        write_reshaped_nifty(rsa_vec_2, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_RSA_%s',ss,anat_names{ii}));
        write_reshaped_nifty(rsa_vec_1-rsa_vec_2, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_peak-RSA_%s',ss,anat_names{ii}));
        write_reshaped_nifty(rsa_vec_2-rsa_vec_1, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_RSA-peak_%s',ss,anat_names{ii}));
    end
    close(h)
end

% rsa_P1 = cellfun(@(x) (sum(x> r2t(0.001,length(A1_corr_vec)))./length(x))*100,rsa_P1_);
% 
% S_mat = squeeze(mean(rsa_P1));
% S_err = squeeze(std(rsa_P1))./sqrt(3);
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
% yline(r2t(0.05,length(M1_anat_vec)))
% xticks(1:nanat)
% xticklabels(anat_names);
% % legend({'Perceptual','Chemical','Mutual'})
% % legend()
% % Subject data points
% c_s = {'r','g','b'}; % Data dots for subjects
% for ii = 1:nanat % For bars for perceptual, chemical and combinations
%     for jj = 1:3
%         plot(x_m(:,ii),squeeze(rsa_P1(jj,ii,:)),c_s{jj},'handle','off')
%     end
% end
% r2t(0.05,length(A2_corr_vec))
% legend('S_I','S_O','sniff col.','neural collinearity')
% ylabel('Representational Similarity (%)')
% % yline(r2t(0.05,sum(utl_mask2(:))));
% % yline(r2t(0.05,nchoosek(length( group_vec),2)));
% savefig(fullfile(savepath,'feat_map'))
% print(fullfile(savepath,'feat_map'),'-dpng')

% clear Fmat_1_m behav_mat unity M_anat M_reg n_reg M_unstack fless_mat fless_mat_unn modelmd_ modelmd S_omat_vals utl_mask utl_mask2
save(fullfile(savepath,'ARC_RSA'),'settings_','rsa_vec_1','rsa_vec_2')

