 %% General Settings
root = 'C:\Work';
tic
settings_.nodor = 160;
settings_.wind = 7500; % Number of samples

% m_id = [8 13; 8 15; 11 16]; % N-settings_.wind
% m_id = [13 61; 15 42; 16 41]; % N-settings_.wind
% m_id = [3 4; 3 4; 3 4]; % INhale Exhale

settings_.nsniffcomp = 14; % 14 for minimal space, 31 for full space
settings_.featorflessnot = true;
settings_.pcamaker = false;

dirs = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
    fullfile(root ,'\SFP\sfp_behav_s02_correct');
    fullfile(root ,'\SFP\sfp_behav_s04_correct')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
    fullfile(root,'ARC\ARC\ARC02\single');
    fullfile(root,'ARC\ARC\ARC03\single')};

savepath = 'C:\Work\SFP\Clustering\Feat_main_updated_fless';
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
rsa_P1 = zeros(3,nanat,length(settings_.loadvec),2); % Complete
hold on
anat_cell = {};
% Subject - index
behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));

hold on
sbplt = 0;
coeffmat = cell(3,1);
figure()
n_vox = zeros(3,nanat,2);
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


    if settings_.pcamaker
        [coeff,sc,~,~,var] = pca(mainmat);
        coeffmat{ss} = coeff;
        
        varsum = cumsum(var);
        varmat{ss} = varsum;
        vardiff{ss} = var;
        subplot(1,3,ss)
        hold on
        plot(varsum)
        xlabel('Num descriptors')
        ylabel('PCA variance explained')
    end


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

        utl_mask = logical(triu(ones(length(group_vec)),1)); % All possible odors

        raw_sniff_cont = spm_read_vols(spm_vol( sprintf('SFP%02d_sniff.nii',ss)));
        raw_sniff_cont = raw_sniff_cont(logical(anatmasks(:,:,:,ii)));
        assert(length(raw_sniff_cont)==size(S_omat_vals_r,1));
        thr = tinv(0.95,nchoosek(length(group_vec),2));
        thr = 0;

        mod_sniff_cont = spm_read_vols(spm_vol( sprintf('SFP%02d_perc.nii',ss)));
        mod_sniff_cont = mod_sniff_cont(logical(anatmasks(:,:,:,ii)));
        assert(length(mod_sniff_cont)==size(S_omat_vals_r,1));
        % if and(ss==2,ii==2)          
        %     error('beepboop')
        % end

        raw_vox = and(raw_sniff_cont>thr, mod_sniff_cont<thr);
        mod_vox = and(raw_sniff_cont>thr, mod_sniff_cont>thr);
        n_vox(ss,ii,:) = [sum(mod_vox) sum(raw_vox)];
        if min(sum(raw_vox),sum(mod_vox))>4
            min(sum(raw_vox),sum(mod_vox))
            A2_corr = SFP_extractUpperTriangles(mainmat,utl_mask);
            
            M2_anat = corrcoef(S_omat_vals_r( raw_vox,:));
            [wt,t_sc] = ARC_multicomputeWeights_tsc([A2_corr task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)], M2_anat(utl_mask));
            rsa_P1(ss,ii,:,1) = wt(2:end-3);

            M2_anat = corrcoef(S_omat_vals_r( mod_vox,:));
            [wt,t_sc] = ARC_multicomputeWeights_tsc([A2_corr  task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)], M2_anat(utl_mask));
            rsa_P1(ss,ii,:,2) = wt(2:end-3);
           
        end
    end
end

%% Figure
rsa_pvals = tcdf(rsa_P1,sum(utl_mask,'all'),'upper');

load('C:\Work\SFP\common_datafiles\snifflabels.mat');
num_desc = 8;
narea = 3;
figure()
hold on
rsa_P1(2,1,:,1) = nan;
for nn = 1:narea
    if nn==1
        ylabel('Beta weight')
    end
    subplot(1,narea,nn)
    hold on
    ARC_barplot(abs(squeeze(rsa_P1(:,nn,1:num_desc,:))),false)
    xticks(1:num_desc)
    xticklabels(proper_list(1:num_desc))
    xtickangle(90)
    title(sprintf('Anat: %s',anat_names{nn}))
end

S_mat = squeeze(mean(rsa_P1));
rsa_pvals_mean = tcdf(S_mat,sum(utl_mask,'all'),'upper');

savefig(fullfile(savepath,'feat_weights'))
print(fullfile(savepath,'feat_weights'),'-dpng')
clearvars fless_mat mainmat Fless_mat_pruned Fless_mat unity task_run set_run sess_run anat_cell

% clear Fmat_1_m behav_mat unity M_anat M_reg n_reg M_unstack fless_mat fless_mat_unn modelmd_ modelmd S_omat_vals utl_mask utl_mask2
save(fullfile(savepath,'SFP_snifffeat'),'settings_','rsa_P1')

%% Main fig
narea = 5;
ndisc = [1 5 7];
figure()
hold on
kk = 0;
for nn = ndisc 
    kk = kk+1;
    if nn==1
        ylabel('Beta weight')
    end
    subplot(1,length(ndisc),kk)
    hold on
    ARC_barplot(abs(squeeze(rsa_P1(:,:,nn,:))),false)
    xticks(1:nanat)
    xticklabels(anat_names)
    xtickangle(90)
    ylim([0 inf])
end

S_mat = squeeze(mean(rsa_P1));
rsa_pvals_mean = tcdf(S_mat,sum(utl_mask,'all'),'upper');

savefig(fullfile(savepath,'feat_weights2'))
print(fullfile(savepath,'feat_weights2'),'-dpng')


%% Principal components
% coeffmat2 = cell(3,1); for ss=1:3; coeffmat2{ss} = coeffmat{ss}(:,1:7); end
% 
% [coeffmat_sort,sortod] = SFP_GaleShapley(coeffmat2);
% vardiffsort = cell(3,1); for ss = 2:3; vardiffsort{ss}= vardiff{ss}(sortod{ss-1}); end
% vardiffsort{1} = vardiff{1};
% num_desc = 25;
% figure()
% hold on
% kk = 0;
% for ss = 1:3
%     for pp = 1:7
%         kk = kk+1;
%         subplot(3,7,kk)
%         bar(coeffmat_sort{ss}(:,pp))
%        title(sprintf('Var explained: %.2f',vardiffsort{ss}(pp)))
%         if ss==3
%             xticks(1:2:num_desc)
%             xticklabels(proper_list(1:2:num_desc))
%             xtickangle(90)
%         end
%     end
% end

%% Mean weights
meanw = squeeze(mean(abs(rsa_P1),3,'omitnan'));
ARC_barplot(meanw)
xticks(1:5)
xticklabels(anat_names)
legend({'Unmod','Mod'})
ylabel('Avg. beta weights')
print(fullfile(savepath,'avg_feat_weights'),'-dpng')
% xtickangle(90)