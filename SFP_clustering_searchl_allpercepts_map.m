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
settings_.featorflessnot = true;
settings_.numClusters = 60;
settings_.pcamaker = true;
settings_.mapper = true;

dirs = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
    fullfile(root ,'\SFP\sfp_behav_s02_correct');
    fullfile(root ,'\SFP\sfp_behav_s04_correct')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
    fullfile(root,'ARC\ARC\ARC02\single');
    fullfile(root,'ARC\ARC\ARC03\single')};

savepath = 'C:\Work\SFP\Clustering\perceptual_clust2';
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3.nii';

anat_names = {'PC','AMY','OFC','OT','AON'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwOT.nii','rwAON.nii'};
nanat = length(anat_names);

settings_.fmasker  = true; % Turn on the functional mask
settings_.single_n = false; % Noisepool
settings_.single_c = true; % Cutoff from sign voxels
settings_.mapper   = true;
settings_.loadvec  = [3 4 9:settings_.nsniffcomp ];

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
    tnvox = sum(anatmasks,'all')*18;

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
     kvox = 0;
    for pp = 1:18

        [idx_int,~] = kmeans(behav_ratings(:,pp),settings_.numClusters,'MaxIter',1000,'Replicates',1000);
        Aint = SFP_splitapply_mean(mainmat,idx_int);
        Aint_corr = corrcoef(Aint');

        % Nuisance regressors
        task_run1 = SFP_splitapply_mean(SFP_splitapply_mean(task_run',idx_int)',idx_int);
        sess_run1 = SFP_splitapply_mean(SFP_splitapply_mean(sess_run',idx_int)',idx_int);
        set_run1 = SFP_splitapply_mean(SFP_splitapply_mean(set_run',idx_int)',idx_int);

        %% Representational connectivity
       
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

            rsa_vec_1 = zeros(nvox,1);
            rsa_vec_2 = zeros(nvox,1);
            rsa_vec_3 = zeros(nvox,1);

            for cnt2 = 1:nvox
                indexindex2 = radius_index + lin_index(cnt2);
                indexindex2 = intersect(lin_index, indexindex2);
                indexindex2 = linindexconv(indexindex2);
                S_omat_vals_r = S_omat_vals(indexindex2,:);
                [r1,~] = find(isnan(S_omat_vals_r));
                S_omat_vals_r(r1,:) = [];

                if (size(S_omat_vals_r,1)>1)

                    M1 = SFP_splitapply_mean(S_omat_vals_r',idx_int);
                    M1 = zscore(M1,[],1);
                  

                    M1_anat = corrcoef(M1');
                

                    utl_mask = logical(triu(ones(length(unique(idx_int))),1)); % All possible odors
                    [~,t_sc1] = ARC_multicomputeWeights_tsc([Aint_corr(utl_mask) task_run1(utl_mask) sess_run1(utl_mask) ], M1_anat(utl_mask));
                    rsa_vec_1(cnt2) = t_sc1(2);

                end
                kvox = kvox+1;
                % waitbar(kvox/(tnvox), h, sprintf('Processing S:%02d, %.2f %%', ss, (kvox/(tnvox)*100))) % Update progress
            end

            rsa_P1_{ss,ii,pp} = rsa_vec_1;
        
            map_area(:,:,:,ii,1)  = unmasker(rsa_vec_1,logical(anatmasks(:,:,:,ii)));

        end

        map_area(map_area==0)=nan;
        map_mat = squeeze(mean(map_area,4,"omitnan"));

        pname = behav.behav(ss).percepts{pp};
        write_reshaped_nifty(squeeze(map_mat(:,:,:,1)), savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_%s',ss,pname));
        % close(h)
    end
    close(h)
end
% rsa_P1 = cellfun(@(x) (sum(x> tinv(0.99,sum(utl_mask,'all')))./length(x))*100,rsa_P1_);
rsa_P1_int = SFP_reshapeCellArray(rsa_P1_);
rsa_P1 = cellfun(@(x) mean(sum(x>1.65,2)),rsa_P1_int);

rsa_P12 = cellfun(@(x) sum(x>1.65,2),rsa_P1_int,'UniformOutput',false);
SFP_performRepeatedMeasuresANOVA(rsa_P12)
% rsa_P1 = cellfun(@(x) median(sum(x>3,2)),rsa_P1_int);

% rsa_P1 =  cellfun(@(x) mean(x),rsa_P1_);
% rsa_P1 = sum(rsa_P1>3,3);
% rsa_P1 =  cellfun(@(x) sum(x>2),rsa_P1_(:,:,1));
% % rsa_P1(3,5,:) = nan;
S_mat = squeeze(nanmean(rsa_P1));
S_err = squeeze(nanstd(rsa_P1))./sqrt(3);
bar(S_mat)
hold on
errorbar([1:5],S_mat,S_err,'.')
xticks([1:5])
xticklabels(anat_names)
for jj = 1:3
    plot(1:5,squeeze(rsa_P1(jj,:)'),'handle','off')
end
anova1(rsa_P1)
save(fullfile(savepath,'ARC_RSA'))
% gcf
% % yline(r2t(0.05,length(M1_anat_vec)))
% xticks(1:nanat)
% xticklabels(anat_names);
% % r2t(0.05,length(A2_corr_vec))
% legend('S_I','S_O')
% ylabel('Representational Similarity (%)')
% savefig(fullfile(savepath,'feat_map'))
% print(fullfile(savepath,'feat_map'),'-dpng')
% 
% rsa_P1_comp = SFP_computePercentages(rsa_P1_,tinv(0.95,sum(utl_mask,'all')));
% ARC_barplot(rsa_P1_comp)
% gcf
% % yline(r2t(0.05,length(M1_anat_vec)))
% xticks(1:nanat)
% xticklabels(anat_names);
% % r2t(0.05,length(A2_corr_vec))
% legend('S_I>S_O','S_O>S_I','Unclassified')
% ylabel('Representational Similarity (%)')
% savefig(fullfile(savepath,'feat_map_comp'))
% print(fullfile(savepath,'feat_map_comp'),'-dpng')
% 
% clearvars fless_mat mainmat Fless_mat_pruned Fless_mat unity task_run set_run sess_run anat_cell
% 
% % clear Fmat_1_m behav_mat unity M_anat M_reg n_reg M_unstack fless_mat fless_mat_unn modelmd_ modelmd S_omat_vals utl_mask utl_mask2
% save(fullfile(savepath,'ARC_RSA'))
% toc
% 
%% Figure Normalization
% Normalize the images into MNI space
img_nrmr = true;
tic
for s = 1:3
    matlabbatch = [];
    dirs = {'\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Neurology\Kahnt_Lab\Vivek\KData\NEMO\NEMO_01\imaging\nii\master_anat\sNEMO01.nii'
        '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Neurology\Kahnt_Lab\Vivek\KData\NEMO\NEMO_02\imaging\nii\anat\sNEMO02.nii'
        '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Neurology\Kahnt_Lab\Vivek\KData\NEMO\NEMO_04\imaging\nii\anat\sNEMO04.nii'};
    wdir = pwd;
%   f2 = dir(fullfile(wdir,sprintf('S%01d',s),'spmT_0006.nii'));
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

%% Image averager
img_avg = true;

if img_avg
    outpurdir = pwd;
    files = SFP_organizeFiles(pwd);

    for ff = 1:size(files,2)
        iminput = files(:,ff);

        tokens = regexp(iminput{1}, 'wSFP(\d+)_(.*).nii', 'tokens');
        strings = tokens{1}{2};  % Extract string part
        fname = strjoin({'gSFP',strings},'_');
        matlabbatch = [];

        matlabbatch{1}.spm.util.imcalc.input = iminput;
        matlabbatch{1}.spm.util.imcalc.output = fname;
        matlabbatch{1}.spm.util.imcalc.outdir = {outpurdir};
%       matlabbatch{1}.spm.util.imcalc.expression = '(double(i1>0.00053)+double(i2>0.00053)+double(i3>0.00053))/3';
%       matlabbatch{1}.spm.util.imcalc.expression = '(-i1-i2-i3)/3';
        % matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3)/3';
        matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run', matlabbatch);
    end
end

