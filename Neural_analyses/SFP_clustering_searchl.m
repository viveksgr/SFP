%% General Settings
root = 'C:\Work';

% Settings
settings_.nodor = 160;
settings_.wind = 4000; % Number of samples
settings_.ncomps = [5 5 5 5];
settings_.nsniffcomp = 31;
settings_.featorflessnot = false;
settings_.numClusters = 60;
settings_.numClusters2 = 60;
settings_.clustbehav = true;
settings_.sniffupdate = false;
settings_.pcamaker = false;
settings_.pcathr = 0.99;
settings_.mapper = true;
settings_.orthog = false;
settings_.depercept = false;

dirs = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
    fullfile(root ,'\SFP\sfp_behav_s02_correct');
    fullfile(root ,'\SFP\sfp_behav_s04_correct')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
    fullfile(root,'ARC\ARC\ARC02\single');
    fullfile(root,'ARC\ARC\ARC03\single')};

savepath = 'C:\Work\SFP\Clustering\Fless_main_updated_orthog';
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
settings_.intplsreg = false;

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')
rsa_P1_ = cell(3,nanat,3); % Complete
rsa_P1t_ = cell(3,nanat,3); % Complete

hold on
anat_cell = {};
% Subject - index
behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));

hold on
sbplt = 0;
numpcs = [13 11 11]; % 90% Variance
idxmats = cell(3,3);
observed_ari= zeros(3,1);
p_value_ari = zeros(3,1);
ari = zeros(3);
thresh = zeros(3,1);

thr_fdr = zeros(3,2);
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

        Fless_mat_pruned = SFP_temporalPCA(Fless_mat_pruned,settings_.pcathr);
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
    dister = 'sqeuclidean';
    if settings_.featorflessnot; mainmat = Feat_mat_pruned; nrep = 1000;  else; mainmat = Fless_mat_pruned; nrep = 10; if ~settings_.pcamaker; dister = 'correlation'; end;  end
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
        [idx,A1] = kmeans(mainmat,settings_.numClusters,'Distance',dister,'MaxIter',1000,'Replicates',nrep);
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
    if settings_.intplsreg
        behav_ratings = behav_ratings(group_vec,3:end);
    else
        behav_ratings = behav_ratings(group_vec,:);
    end


    [idx2,~] = kmeans(behav_ratings,settings_.numClusters2,'MaxIter',1000,'Replicates',nrep);
    if settings_.clustbehav; group_vec = idx2; end
    A2 = SFP_splitapply_mean(mainmat,group_vec);
    A2_corr = corrcoef(A2');

    idx3 = datasample(1:settings_.numClusters,length(idx));
    idxmats{ss,1}= idx;
    idxmats{ss,2}= idx2;
    idxmats{ss,3}= idx3;

    A3 = SFP_splitapply_mean(mainmat,idx3);
    A3_corr = corrcoef(A3');

    [observed_ari(ss),p_value_ari(ss),threshold]=aritester(idx,idx2);
    thresh(ss)= threshold;
    if settings_.orthog
        [idx, idx2,ari(:,ss)] = sfp_decoupleClusteringIndices_prob(idx, idx2, A1_corr, A2_corr, 10000, threshold,3);
    end

    % Nuisance regressors
    task_run1 = SFP_splitapply_mean(SFP_splitapply_mean(task_run',idx)',idx);
    sess_run1 = SFP_splitapply_mean(SFP_splitapply_mean(sess_run',idx)',idx);
    set_run1 = SFP_splitapply_mean(SFP_splitapply_mean(set_run',idx)',idx);
    task_run2 = SFP_splitapply_mean(SFP_splitapply_mean(task_run',group_vec)',group_vec);
    sess_run2 = SFP_splitapply_mean(SFP_splitapply_mean(sess_run',group_vec)',group_vec);
    set_run2 = SFP_splitapply_mean(SFP_splitapply_mean(set_run',group_vec)',group_vec);
    task_run3 = SFP_splitapply_mean(SFP_splitapply_mean(task_run',idx3)',idx3);
    sess_run3 = SFP_splitapply_mean(SFP_splitapply_mean(sess_run',idx3)',idx3);
    set_run3 = SFP_splitapply_mean(SFP_splitapply_mean(set_run',idx3)',idx3);

    if  settings_.depercept
        b1 = SFP_splitapply_mean(behav_ratings,idx);
        b1corr = corrcoef(b1');
        b2 = SFP_splitapply_mean(behav_ratings,group_vec);
        b2corr = corrcoef(b2');
        u1corr = SFP_splitapply_mean(SFP_splitapply_mean(unity',idx)',idx);
        u2corr = SFP_splitapply_mean(SFP_splitapply_mean(unity',group_vec)',group_vec);
    end
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

        rsa_vec_1w = zeros(nvox,1);
        rsa_vec_2w = zeros(nvox,1);
        rsa_vec_3w = zeros(nvox,1);

        rsa_vec_1t = zeros(nvox,1);
        rsa_vec_2t = zeros(nvox,1);
        rsa_vec_3t = zeros(nvox,1);

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
                M3 = SFP_splitapply_mean(S_omat_vals_r',idx3);
                M3 = zscore(M3,[],1);
                M1_anat = corrcoef(M1');
                M2_anat = corrcoef(M2');
                M3_anat = corrcoef(M3');

                if ~settings_.depercept
                    utl_mask = logical(triu(ones(length(unique(idx))),1)); % All possible odors
                    [wt1,t_sc1] = ARC_multicomputeWeights_tsc([A1_corr(utl_mask) task_run1(utl_mask) sess_run1(utl_mask) set_run1(utl_mask)], M1_anat(utl_mask));
                    rsa_vec_1w(cnt2) = wt1(2);
                    rsa_vec_1t(cnt2) = t_sc1(2);

                    utl_mask = logical(triu(ones(length(unique(group_vec))),1)); % All possible odors
                    [wt2,t_sc2] = ARC_multicomputeWeights_tsc([A2_corr(utl_mask) task_run2(utl_mask) sess_run2(utl_mask)], M2_anat(utl_mask));
                    rsa_vec_2w(cnt2) = wt2(2);
                    rsa_vec_2t(cnt2) = t_sc2(2);

                    utl_mask = logical(triu(ones(length(unique(idx3))),1)); % All possible odors
                    [wt3,t_sc3] = ARC_multicomputeWeights_tsc([A3_corr(utl_mask) task_run3(utl_mask) sess_run3(utl_mask) set_run3(utl_mask)], M3_anat(utl_mask));
                    rsa_vec_3w(cnt2) = wt3(2);
                    rsa_vec_3t(cnt2) = t_sc3(2);
                else
                    utl_mask = logical(triu(ones(length(unique(idx))),1)); % All possible odors

                    [wt1,t_sc1] = ARC_multicomputeWeights_tsc([A1_corr(utl_mask) b1corr(utl_mask) u1corr(utl_mask)  task_run1(utl_mask) sess_run1(utl_mask) set_run1(utl_mask)], M1_anat(utl_mask));
                    % [~,t_sc] = ARC_multicomputeWeights_tsc([A1_corr(utl_mask)], M1_anat(utl_mask));

                    rsa_vec_1(cnt2) = wt1(2);

                    utl_mask = logical(triu(ones(length(unique(group_vec))),1)); % All possible odors
                    % set_run2 and sess_run2 are collinear
                    [wt2,t_sc2] = ARC_multicomputeWeights_tsc([A2_corr(utl_mask) b2corr(utl_mask) task_run2(utl_mask) sess_run2(utl_mask)], M2_anat(utl_mask));
                    % [~,t_sc] = ARC_multicomputeWeights_tsc([A2_corr(utl_mask)], M2_anat(utl_mask));

                    rsa_vec_2(cnt2) = wt2(2);
                end

            end
            kvox = kvox+1;
            waitbar(kvox/(tnvox), h, sprintf('Processing S:%02d, %.2f %%', ss, (kvox/(tnvox)*100))) % Update progress
        end


        if and(ii==3,ss==2)
            'Arre bitwa...'
        end

        rsa_P1_{ss,ii,1} = rsa_vec_1w;
        rsa_P1_{ss,ii,2} = rsa_vec_2w;
        rsa_P1_{ss,ii,3} = rsa_vec_3w;
        rsa_P1t_{ss,ii,1} = rsa_vec_1t;
        rsa_P1t_{ss,ii,2} = rsa_vec_2t;
        rsa_P1t_{ss,ii,3} = rsa_vec_3t;

        map_area(:,:,:,ii,1) = unmasker(rsa_vec_1t,logical(anatmasks(:,:,:,ii)));
        map_area(:,:,:,ii,2) = unmasker(rsa_vec_2t,logical(anatmasks(:,:,:,ii)));
        map_area(:,:,:,ii,3) = unmasker(rsa_vec_3t,logical(anatmasks(:,:,:,ii)));
        %
        % if settings_.mapper
        %     rsa_vec_1 = unmasker(rsa_vec_1,logical(anatmasks(:,:,:,ii)));
        %     rsa_vec_2 = unmasker(rsa_vec_2,logical(anatmasks(:,:,:,ii)));
        %     write_reshaped_nifty(rsa_vec_1, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_peak_%s',ss,anat_names{ii}));
        %     write_reshaped_nifty(rsa_vec_2, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_RSA_%s',ss,anat_names{ii}));
        %     write_reshaped_nifty(rsa_vec_1-rsa_vec_2, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_peak-RSA_%s',ss,anat_names{ii}));
        %     write_reshaped_nifty(rsa_vec_2-rsa_vec_1, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_RSA-peak_%s',ss,anat_names{ii}));
        % end

    end

    map_area(map_area==0)=nan;
    map_mat = squeeze(mean(map_area,4,"omitnan"));

    m1 = squeeze(map_mat(:,:,:,1));
    m2 = squeeze(map_mat(:,:,:,2));
    m3 = squeeze(map_mat(:,:,:,3));

    df1 = sum(~isnan(m1(:)));
    func = @(x) 2 * (1 - tcdf(abs(x),df1));
    p1_mat = arrayfun(func,m1(~isnan(m1)));

    df2 = sum(~isnan(m2(:)));
    func = @(x) 2 * (1 - tcdf(abs(x),df2));
    p2_mat = arrayfun(func,m2(~isnan(m2)));

    df3 = sum(~isnan(m3(:)));
    func = @(x) 2 * (1 - tcdf(abs(x),df3));
    p3_mat = arrayfun(func,m3(~isnan(m3)));


    thr_fdr(ss,1) = tinv(1-fdr_benjhoc(p1_mat),df1);
    thr_fdr(ss,2) = tinv(1-fdr_benjhoc(p2_mat),df2);
    thr_fdr(ss,3) = tinv(1-fdr_benjhoc(p3_mat),df3);


    write_reshaped_nifty(squeeze(map_mat(:,:,:,1)), savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_sniff',ss));
    write_reshaped_nifty(squeeze(map_mat(:,:,:,2)), savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_perc',ss));
    write_reshaped_nifty(squeeze(map_mat(:,:,:,3)), savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_rand',ss));
    close(h)
end

% rsa_P1 = cellfun(@(x) (sum(x> tinv(0.95,sum(utl_mask,'all')))./length(x))*100,rsa_P1_);
rsa_P1 =  cellfun(@(x) mean(x),rsa_P1_);
rsa_P1t = cellfun(@(b, t) SFP_calculateAvgTscore(b, t), rsa_P1_, rsa_P1t_, 'UniformOutput', true);
df = sum( utl_mask,'all');


[adjustedPValues, t_st,df] = SFP_performRepeatedMeasuresANOVA_diff(rsa_P1_)


[diffMat1, diffMat2, diffMat3] = SFP_transformCellArray(rsa_P1_);
pValue = SFP_performRepeatedMeasuresANOVA(diffMat3)
p_values_matrix = ARC_RSA_pvals(rsa_P1t, rsa_P1, df)

% rsa_P1(3,5,:) = nan;
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

SFP_clearLargeVariables
save(fullfile(savepath,'ARC_RSA'))
toc

%% Figure Normalization
% Normalize the images into MNI space
cd(savepath)
img_nrmr = true;
tic
for s = 1:3
    matlabbatch = [];
    dirs = {'\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Neurology\Kahnt_Lab\Vivek\KData\NEMO\NEMO_01\imaging\nii\master_anat\sNEMO01.nii'
        '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Neurology\Kahnt_Lab\Vivek\KData\NEMO\NEMO_02\imaging\nii\anat\sNEMO02.nii'
        '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Neurology\Kahnt_Lab\Vivek\KData\NEMO\NEMO_04\imaging\nii\anat\sNEMO04.nii'};
    wdir = pwd;
    %     f2 = dir(fullfile(wdir,sprintf('S%01d',s),'spmT_0006.nii'));
    f2 = dir(fullfile(wdir,sprintf('SFP%02d*.nii',s)));

    files = {};
    for zz = 1:length(f2); files{zz,1} = fullfile(wdir,f2(zz).name); end

    if img_nrmr
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {dirs{s}};
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = files;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'C:\Work\Tools\spm12\tpm\TPM.nii'};
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
        spm_jobman('run', matlabbatch);
    end
end
toc

%% Image averager
img_avg = true;

if img_avg
    outpurdir = pwd;
    files = SFP_organizeFiles(pwd);

    for ff = 1:size(files,2)
        iminput = files(:,ff);

        tokens = regexp(iminput{1}, 'wSFP(\d+)_(.*).nii', 'tokens');
        strings = tokens{1}{2};  % Extract string part
        fname = strjoin({'cgSFP',strings},'_');
        matlabbatch = [];

        matlabbatch{1}.spm.util.imcalc.input = iminput;
        matlabbatch{1}.spm.util.imcalc.output = fname;
        matlabbatch{1}.spm.util.imcalc.outdir = {outpurdir};
        if ff==1
            matlabbatch{1}.spm.util.imcalc.expression = '(double(i1>1.79)+double(i2>1.65)+double(i3>2.37))';
        else
            matlabbatch{1}.spm.util.imcalc.expression = '(double(i1>2.34)+double(i2>1.70)+double(i3>2.09))';
        end
        % matlabbatch{1}.spm.util.imcalc.expression = '(-i1-i2-i3)/3';
        % matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3)/3';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run', matlabbatch);
    end
end

