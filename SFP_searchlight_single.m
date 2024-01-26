%% General Settings
root = 'C:\Work';
nodor = 160;
wind = 75; % Number of samples

m_id = [8 13; 8 15; 11 16]; % N-wind
% m_id = [13 61; 15 42; 16 41]; % N-wind
% m_id = [3 4; 3 4; 3 4]; % INhale Exhale

dirs = {fullfile(root ,'\SFP\sfp_behav_s01');
        fullfile(root ,'\SFP\sfp_behav_s02');
        fullfile(root ,'\SFP\sfp_behav_s03')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
        fullfile(root,'ARC\ARC\ARC02\single');
        fullfile(root,'ARC\ARC\ARC03\single')};

dirs3 = {fullfile(root,'SFP\basic_peaks\SFP1');
        fullfile(root,'SFP\basic_peaks\SFP2');
        fullfile(root,'SFP\basic_peaks\SFP3')};

maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3.nii'; 

anat_names = {'PC','AMY','OFC','OT','AON'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwOT.nii','rwAON.nii'};
nanat = length(anat_names);

fmasker = true; % Turn on the functional mask
single_n = false; % Noisepool
single_c = true; % Cutoff from sign voxels
mapper = true;

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')
rsa_P1_ = cell(3,nanat,2); % Complete
hold on
anat_cell = {};
% Subject - index

for ss = [1 2 3] % Subject
    fprintf('Subject: %02d\n',ss)
    if ss==3; s2 = 4; else; s2 = ss; end
    statpath = dirs{ss};
    anatdir = dirs2{ss};
    savepath = dirs3{ss};
    mkdir(savepath)


    load(fullfile(statpath,'sfp_feats.mat'))

    % Fless_mat_pruned = vertcat(feat_mat{:});
    % Fmat_1 = Fless_mat_pruned(:,m_id(ss,1));
    % Fmat_2 = Fless_mat_pruned(:,m_id(ss,2));

    Fless_mat = vertcat(fless_mat{:});
    Fless_mat_pruned = Fless_mat(:,1:wind);
    Fmat_1 = Fless_mat_pruned(:,m_id(ss,1));
    Fmat_2 = Fless_mat_pruned(:,m_id(ss,2));


    Fmat_1_m = -abs(Fmat_1-Fmat_1');
    Fmat_2_m = -abs(Fmat_2-Fmat_2');

    %% Gray Matter, Functional and Anatomical Masks
    mask = (spm_read_vols(spm_vol(fullfile(anatdir, maskfile)))); % Mask used to construct odor files
    mask(isnan(mask))=0;
    mask = logical(mask);
    fmask = (spm_read_vols(spm_vol(fullfile(anatdir, fmaskfile)))); % Mask used to examine voxels in RSA
    fmask(isnan(fmask))=0;
    if ~fmasker
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

    %% Defining RDMs
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
    
    %% Representational connectivity
    for ii = 1:length(anat_names)
        fprintf('area:%02d\n',ii)
        modelmd_ = load(fullfile(anatdir,'desniff',anat_names{ii},'TYPEC_FITHRF_GLMDENOISE.mat'),'modelmd','noisepool');
        modelmd = squeeze(modelmd_.modelmd);
        noisepool = modelmd_.noisepool;
        if single_c
            modelmd = modelmd(masks_set_cell{ii},:);
            noisepool = noisepool(masks_set_cell{ii});
        end
        S_omat_vals = zscore(modelmd,[],2);
        [r1,~] = find(isnan(S_omat_vals));
        S_omat_vals(r1,:) = [];
        %         M_anat = 1-corrcoef(S1_omat_vals);
        
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
                temp = corrcoef(M_reg,Fmat_1_m(utl_mask));
                rsa_vec_1(cnt2) = temp(2);
                temp = corrcoef(M_reg,Fmat_2_m(utl_mask));
                rsa_vec_2(cnt2) = temp(2);
            end
        end
        rsa_P1_{ss,ii,1} = rsa_vec_1;
        rsa_P1_{ss,ii,2} = rsa_vec_2;
        
        if mapper
            rsa_vec_1 = unmasker(rsa_vec_1,logical(anatmasks(:,:,:,ii)));
            rsa_vec_2 = unmasker(rsa_vec_2,logical(anatmasks(:,:,:,ii)));
            write_reshaped_nifty(rsa_vec_1, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_peak_%s',ss,anat_names{ii}));
            write_reshaped_nifty(rsa_vec_2, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_RSA_%s',ss,anat_names{ii}));
            write_reshaped_nifty(rsa_vec_1-rsa_vec_2, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_peak-RSA_%s',ss,anat_names{ii}));
            write_reshaped_nifty(rsa_vec_2-rsa_vec_1, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_RSA-peak_%s',ss,anat_names{ii}));
        end
    end
end
clear Fmat_1_m Fmat_2_m unity M_anat M_reg n_reg M_unstack fless_mat fless_mat_unn modelmd_ modelmd S_omat_vals utl_mask utl_mask2
save(fullfile(savepath,'ARC_RSA'))

% rsa_P1_pvals = cellfun(@(x) r2p(x,nchoosek(ntrials,2)),rsa_P1_,'UniformOutput',false);
% rsa_P1_pvals2 = cat(1,rsa_P1_pvals{:});
% rsa_P1_thresh = fdr_benjhoc(rsa_P1_pvals2);

% for pp = 1:numel(rsa_P1_thresh);  if isempty(rsa_P1_thresh{pp}); rsa_P1_thresh{pp} = 0; end; end
% rsa_P1_thresh{1,1,1} = 0; % Check this manually
% rsa_P1_thresh{3,1,2} = 0; % Check this manually

% rsa_P1 = cellfun(@(x) (sum(x>r2t(0.05,nchoosek(ntrials,2)))./length(x))*100,rsa_P1_);
% rsa_P1 = cellfun(@(x1,x2) (sum(x1<x2.*ones(size(x1)))./length(x1))*100,rsa_P1_pvals,rsa_P1_thresh);
% rsa_P1 = cellfun(@(x) (sum(x<rsa_P1_thresh)./length(x))*100,rsa_P1_pvals);
rsa_P1 = cellfun(@(x) mean(x),rsa_P1_);
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
legend('peak','rsa')
ylabel('Representational Similarity (r)')
% yline(r2t(0.05,sum(utl_mask2(:))));
yline(r2t(0.05,nchoosek(length( group_vec),2)));
savefig(fullfile(savepath,'ARC_RSA'))
print(fullfile(savepath,'ARC_RSA'),'-dpng')

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
img_nrmr = false;
tic
for s = 1:3
    matlabbatch = [];
    dirs = {'C:\Data\NEMO\NEMO_01\imaging\nii\master_anat\meansNEMO_01_set1_sess1-0015-00001-000176-01.nii'
        'C:\Data\NEMO\NEMO_02\imaging\nii\anat\sNEMO02.nii'
        'C:\Data\NEMO\NEMO_04\imaging\nii\anat\sNEMO04.nii'};
    wdir = pwd;
%     f2 = dir(fullfile(wdir,sprintf('S%01d',s),'spmT_0006.nii'));
    f2 = dir(fullfile(wdir,sprintf('sfp_behav_s%02d',s),'rsa_vec_fless.nii'));

    files = {};
    for zz = 1:length(f2); files{zz,1} = fullfile(wdir,sprintf('sfp_behav_s%02d',s),f2(zz).name); end
    
    if img_nrmr
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {dirs{s}};
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = files;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'C:\Toolboxes\spm12\tpm\TPM.nii'};
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
img_avg = false;
ntrials = 4320;
thr = r2t(0.05,nchoosek(ntrials,2));
if img_avg
    outpurdir = 'C:\Data\SFP\neural';
    dirs = {'C:\Data\SFP\sfp_behav_s01'
            'C:\Data\SFP\sfp_behav_s02'
            'C:\Data\SFP\sfp_behav_s03'};
%           outpurdir = 'C:\Data\ARC\ARC\val_glm\val-sal\basic';
%     dirs = {'C:\Data\ARC\ARC\val_glm\val-sal\basic\S1'
%             'C:\Data\ARC\ARC\val_glm\val-sal\basic\S2'
%             'C:\Data\ARC\ARC\val_glm\val-sal\basic\S3'};
    files = {};
    for ss = 1:3
        f2 = dir(fullfile(dirs{ss},'wrsa_vec_fless.nii'));
        for zz = 1:length(f2); files{zz,ss} = fullfile(dirs{ss},f2(zz).name); end
    end
    
    for ff = 1:size(files,1)
        iminput = files(ff,:)';
        fname = files{ff,1}(end-8:end-4);
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