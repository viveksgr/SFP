%% General Settings
root = 'C:\Work';
settings_.nodor = 160;
settings_.wind = 75; % Number of samples

% m_id = [8 13; 8 15; 11 16]; % N-settings_.wind
% m_id = [13 61; 15 42; 16 41]; % N-settings_.wind
% m_id = [3 4; 3 4; 3 4]; % INhale Exhale
settings_.ncomps = [5 5 5 5];
settings_.nsniffcomp = 14;
settings_.ncomp_plt = 8;
settings_.prfm_pca = false;
settings_.reg_mat = true;

dirs = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
        fullfile(root ,'\SFP\sfp_behav_s02_correct');
        fullfile(root ,'\SFP\sfp_behav_s04_correct')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
        fullfile(root,'ARC\ARC\ARC02\single');
        fullfile(root,'ARC\ARC\ARC03\single')};

dirs3 = {fullfile(root ,'\SFP\sfp_behav_s01_correct_reg');
        fullfile(root ,'\SFP\sfp_behav_s02_correct_reg');
        fullfile(root ,'\SFP\sfp_behav_s04_correct_reg')};


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

figure('Position',[0.5 0.5 1280 360])
hold on

sbplt = 0;
for ss = [1 2 3] % Subject
    h = waitbar(0, 'Please wait...'); % Progress bar per subject
    fprintf('Subject: %02d\n',ss)
    if ss==3; s2 = 4; else; s2 = ss; end
    statpath = dirs{ss};
    anatdir = dirs2{ss};
    savepath = dirs3{ss};
    mkdir(savepath)

    load(fullfile(statpath,'sfp_feats_main.mat'))
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
    [~,argsort] = sort(vertcat(onsets{:}));
    group_vec = group_vec(argsort);
    unity = unity(argsort,argsort);
    utl_mask = logical(triu(ones(length(unity)),1)); % All possible odors

    % Behavioral RSMs
    behav_ratings = behav.behav(ss).ratings;
    behav_ratings = behav_ratings(group_vec,:);
    behav_mat = corrcoef(behav_ratings');
    
    % CCA
    % Assuming F1 is 4320x8 and F2 is 4320x18
    % Remove constant columns from F1
    F1 = Feat_mat_pruned;
    F2 = behav_ratings;
    F1(isnan(F1))=0;
    F2(isnan(F2))=0;

    % Standardize the data
    F1_std = zscore(F1);
    F2_std = zscore(F2);

    % Remove near-constant columns from F1
    varF1 = var(F1_std);
    F1_no_near_constant = F1_std(:, varF1 > 1e-10); % Threshold for near-constant

    % Remove near-constant columns from F2
    varF2 = var(F2_std);
    F2_no_near_constant = F2_std(:, varF2 > 1e-10); % Threshold for near-constant


    % Apply canonical correlation analysis
    [F1_coeff, F2_coeff, r, U,V,stats] = canoncorr(F1_no_near_constant, F2_no_near_constant);


    


    
    for plt = 1:settings_.ncomp_plt
        sbplt = sbplt+1;
        subplot(3,settings_.ncomp_plt,sbplt) 
        bar(F1_coeff(:,plt))    
        % xticklabels({'Inflow','ExFlow','Time Peak','Time Trough','InVol','ExVol','InDur','ExDur'})
        title(sprintf('CCA%01d r:%.2f',plt,r(plt)))
    end

    U = U(:,1:settings_.ncomps(ss));

    if settings_.prfm_pca
        [~,A1_mat] = pca(F1_std);
        A1_mat = A1_mat(:,1:settings_.ncomps(ss));
    else
        if settings_.reg_mat
            A1_mat = SFP_multiregressmeout(F1_std, F2_no_near_constant);
        else
            A1_mat = F1_std;
        end
    end


    A1 = corrcoef(A1_mat');
    A2 = corrcoef(U');

    kvox = 0;
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
        S_omat_vals = zscore(modelmd,[],2);
        [r1,~] = find(isnan(S_omat_vals));
        S_omat_vals(r1,:) = [];
        % M_anat = 1-corrcoef(S1_omat_vals);

        % Searchlight configuration
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
            M_anat = corrcoef(S_omat_vals_r);

            if (size(S_omat_vals_r,1)>1)
                M_unstack = M_anat(logical(utl_mask));
                n_reg = unity(logical(utl_mask));
                M_reg = regressmeout(M_unstack',n_reg')';

                [~,temp] = SFP_computeWeights_tsc( A1(utl_mask),  A2(utl_mask),M_reg);
                rsa_vec_1(cnt2) = temp(2);
                rsa_vec_2(cnt2) = temp(3);
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
savefig(fullfile(savepath,'CCA'))
print(fullfile(savepath,'CCA.png'),'-dpng')


% rsa_P1_pvals = cellfun(@(x) r2p(x,nchoosek(ntrials,2)),rsa_P1_,'UniformOutput',false);
% rsa_P1_pvals2 = cat(1,rsa_P1_pvals{:});
% rsa_P1_thresh = fdr_benjhoc(rsa_P1_pvals2);

% for pp = 1:numel(rsa_P1_thresh);  if isempty(rsa_P1_thresh{pp}); rsa_P1_thresh{pp} = 0; end; end
% rsa_P1_thresh{1,1,1} = 0; % Check this manually
% rsa_P1_thresh{3,1,2} = 0; % Check this manually

% rsa_P1 = cellfun(@(x) (sum(x>r2t(0.05,nchoosek(ntrials,2)))./length(x))*100,rsa_P1_);
% rsa_P1 = cellfun(@(x1,x2) (sum(x1<x2.*ones(size(x1)))./length(x1))*100,rsa_P1_pvals,rsa_P1_thresh);
rsa_P1 = cellfun(@(x) (sum(x> tinv(1 - 0.001, nchoosek(size(Fless_mat,1),2)))./length(x))*100,rsa_P1_);
% rsa_P1 = cellfun(@(x) mean(x),rsa_P1_);


thr = tinv(1 - 0.05, nchoosek(size(Fless_mat,1),2));
rsa_P1 = SFP_computePercentages(rsa_P1_,thr);

S_mat = squeeze(mean(rsa_P1));
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
legend('Sniff','Perceptual')
ylabel('Representational Similarity (%)')
% yline(r2t(0.05,sum(utl_mask2(:))));
% yline(r2t(0.05,nchoosek(length( group_vec),2)));
savefig(fullfile(savepath,'ARC_RSA_sign_comp_neg'))
print(fullfile(savepath,'ARC_RSA_sign_comp_neg'),'-dpng')

% clear Fmat_1_m behav_mat unity M_anat M_reg n_reg M_unstack fless_mat fless_mat_unn modelmd_ modelmd S_omat_vals utl_mask utl_mask2
save(fullfile(savepath,'ARC_RSA'),'settings','rsa_P1') 

%% Bar plots of Valence and Salience as a function of other dims
barplotter = false;

if barplotter
    chemer = false;
    valier = false;
    niter = 10000;
    nC = 20; % NUmber of chemical components

    figure('Position',[0.5 0.5 1280 360])
    hold on
    nS = 3;
    for ss = 1:nS
        subplot(1,3,ss)
        hold on
        Y = behavP.behav(ss).ratings(:,2);
        if ~valier;  Y = abs(Y); end
        if chemer
            X = behavC.behav(ss).ratings(:,1:nC); 
            xp = [1: nC]; % Labels
        else 
            X = behavP.behav(ss).ratings(:,[1 3:end]); 
            xp = behavP.behav(ss).percepts([1 3:end]);
        end

        [beta_m, beta_err, p_values,~,res] = bootstrapRidgeres(X, Y, niter, 0.1);
        plotBarWithSignificance(beta_m, beta_err*1.96, p_values)
        % Note the p_value is approximate since it's not strictly LOOCV
        title(sprintf('S: %02d r: %.2f p: %.3f',ss,res,r2p(res,160)))
        if  ~chemer; xticks(1:17); xtickangle(90); xticklabels(xp); end
        ylabel('Ridge regression weights')
    end
    % savefig('Sal_perceptual2')
    % print('Sal_perceptual2','-dpng')
end



%% Normalize the images into MNI space
img_nrmr = true;
tic
dirs3 = {fullfile(root ,'\SFP\sfp_behav_s01_cluster');
    fullfile(root ,'\SFP\sfp_behav_s02_cluster');
    fullfile(root ,'\SFP\sfp_behav_s04_cluster')};
for s = 1:3
    matlabbatch = [];
    dirs = {'C:\Work\ARC\ARC\ARC01\masks\sNEMO01.nii'
        'C:\Work\ARC\ARC\ARC02\masks\sNEMO02.nii'
        'C:\Work\ARC\ARC\ARC03\masks\sNEMO04.nii'};

    wdir = dirs3{s};
    %     f2 = dir(fullfile(wdir,sprintf('S%01d',s),'spmT_0006.nii'));
    f2 = dir(fullfile(wdir,'SFP*.nii'));

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

anat_names = {'PC','AMY','OFC','OT','AON'}; 

if img_avg
    outpurdir = 'C:\Work\SFP\Clustering';
    dirs = {'C:\Work\SFP\sfp_behav_s01_cluster'
            'C:\Work\SFP\sfp_behav_s02_cluster'
            'C:\Work\SFP\sfp_behav_s04_cluster'};
%           outpurdir = 'C:\Data\ARC\ARC\val_glm\val-sal\basic';
%     dirs = {'C:\Data\ARC\ARC\val_glm\val-sal\basic\S1'
%             'C:\Data\ARC\ARC\val_glm\val-sal\basic\S2'
%             'C:\Data\ARC\ARC\val_glm\val-sal\basic\S3'};


    files = {};
    for ss = 1:3
        f2 = dir(fullfile(dirs{ss},'w*.nii'));
        for zz = 1:length(f2); files{zz,ss} = fullfile(dirs{ss},f2(zz).name); end
    end
    
    for ff = 1:size(files,1)
        iminput = files(ff,:)';
        fname = files{ff,1}(end-15:end-4);
        matlabbatch = [];
        
        matlabbatch{1}.spm.util.imcalc.input = iminput;
        matlabbatch{1}.spm.util.imcalc.output = fname;
        matlabbatch{1}.spm.util.imcalc.outdir = {outpurdir};
%       matlabbatch{1}.spm.util.imcalc.expression = '(double(i1>0.00053)+double(i2>0.00053)+double(i3>0.00053))/3';
%         matlabbatch{1}.spm.util.imcalc.expression = '(-i1-i2-i3)/3';
         matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3)/3';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run', matlabbatch);
    end
end