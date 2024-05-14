%% General Settings
root = 'C:\Work';
tic
settings_.nodor = 160;
settings_.wind = 7500; % Number of samples
delta = 0.25;
timeaxis = 0:delta:7;
timeaxis2 = [timeaxis 7+delta]*1000;


% m_id = [8 13; 8 15; 11 16]; % N-settings_.wind
% m_id = [13 61; 15 42; 16 41]; % N-settings_.wind
% m_id = [3 4; 3 4; 3 4]; % INhale Exhale

settings_.nsniffcomp = 31; % 14 for minimal space, 31 for full space
settings_.featorflessnot = false;
settings_.multiregress = true;
settings_.pcamaker = false;

dirs = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
    fullfile(root ,'\SFP\sfp_behav_s02_correct');
    fullfile(root ,'\SFP\sfp_behav_s04_correct')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
    fullfile(root,'ARC\ARC\ARC02\single');
    fullfile(root,'ARC\ARC\ARC03\single')};

filterpath = 'C:\Work\SFP\Clustering\Feat_main_updated';
savepath = 'C:\Work\SFP\Feat_temporal';
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3.nii';

% anat_names = {'PC','AMY','OFC','OT','AON'};
% anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwOT.nii','rwAON.nii'};
anat_names = {'PC','AMY','OFC'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii'};

nanat = length(anat_names);

settings_.fmasker = true; % Turn on the functional mask
settings_.single_n = false; % Noisepool
settings_.single_c = true; % Cutoff from sign voxels
settings_.mapper = true;
settings_.loadvec = [3 4 9:settings_.nsniffcomp ];

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')
rsa_P1 = zeros(3,nanat,length(timeaxis),2); % Complete
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
    mkdir(savepath)

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
    behav_corr = corrcoef(behav_ratings');
    utl_mask = logical(triu(ones(length(group_vec)),1)); % All possible odors

    %% Representational connectivity
    for ii = 1:length(anat_names)
        fprintf('area:%02d\n',ii)

        raw_sniff_cont = spm_read_vols(spm_vol( fullfile(filterpath,sprintf('SFP%02d_peak_%s.nii',ss,anat_names{ii}))));
        raw_sniff_cont = raw_sniff_cont(logical(anatmasks(:,:,:,ii)));
        thr = tinv(0.999,nchoosek(length(group_vec),2));

        mod_sniff_cont = spm_read_vols(spm_vol( fullfile(filterpath,sprintf('SFP%02d_RSA_%s.nii',ss,anat_names{ii}))));
        mod_sniff_cont = mod_sniff_cont(logical(anatmasks(:,:,:,ii)));

        raw_vox = and(raw_sniff_cont>thr, mod_sniff_cont<thr);
        mod_vox = and(raw_sniff_cont>thr, mod_sniff_cont>thr);
        n_vox(ss,ii,:) = [sum(mod_vox) sum(raw_vox)];


        for tpoint = 1:length(timeaxis)
            if min(sum(raw_vox),sum(mod_vox))>1

                modelmd_ = load(fullfile(anatdir,'desniff',anat_names{ii},'TYPEC_FITHRF_GLMDENOISE.mat'),'modelmd','noisepool');
                modelmd = squeeze(modelmd_.modelmd);
                % noisepool = modelmd_.noisepool;

                if settings_.single_c
                    modelmd = modelmd(masks_set_cell{ii},:);
                    % noisepool = noisepool(masks_set_cell{ii});
                end
                S_omat_vals_r = modelmd;
                [r1,~] = find(isnan(S_omat_vals_r));
                S_omat_vals_r(r1,:) = [];
                assert(and(length(mod_sniff_cont)==size(S_omat_vals_r,1),length(raw_sniff_cont)==size(S_omat_vals_r,1)));
                % M_anat = 1-corrcoef(S1_omat_vals);
                M2_anat = corrcoef(S_omat_vals_r( raw_vox,:));


                A2_corr = pdist(mainmat(:,timeaxis2(tpoint)+1:timeaxis2(tpoint+1)),"correlation");
                A2_corr = 1-squareform(A2_corr);

                [wt,t_sc] = ARC_multicomputeWeights_tsc([A2_corr(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)], M2_anat(utl_mask));
                rsa_P1(ss,ii,tpoint,1) = wt(2);

                M2_anat = corrcoef(S_omat_vals_r( mod_vox,:));
                [wt,t_sc] = ARC_multicomputeWeights_tsc([A2_corr(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)], M2_anat(utl_mask));
                rsa_P1(ss,ii,tpoint,2) = wt(2);
            end
        end
    end
end

load('C:\Work\SFP\Feat_temporal\timecourse.mat')

% Subject-wise plot
figure('Position',[0 0 1480 720])
kk = 0;
for ss = 1:3
for nn = 1:3
    kk = kk+1;
    subplot(3,3,kk)
    hold on
    yyaxis left

    snifft = interp1((0:length(sniff_trace(ss,:))-1)/10,sniff_trace(ss,:),timeaxis);
    plot(timeaxis,snifft)
    ylabel('Sniff Trace')

    yyaxis right
    % plot((1:wind)/10,mean(corrmod,1))
    plot(timeaxis,squeeze(rsa_P1(ss,nn,:,1)))
    hold on
    plot(timeaxis,squeeze(rsa_P1(ss,nn,:,2)))
    m1 = max(max(squeeze(rsa_P1(ss,nn,4:end,2))),max(squeeze(rsa_P1(ss,nn,4:end,1))));
    m2 = min(min(squeeze(rsa_P1(ss,nn,4:end,2))),min(squeeze(rsa_P1(ss,nn,4:end,1))));
    ylim([m2 m1 ])

    ylabel('RSA performance')
    if and(ss==3,nn==3)
        legend({'Sniff','Unmod voxels','Mod voxels'})
    end
    xlabel('time(s)')
    title(sprintf('RSA for sub %02d, %s',ss,anat_names{nn}))    
end
end

% Subject-wise plot
figure('Position',[0 0 1480 720])
kk = 0;
for ss = 1:3
for nn = 1:3
    kk = kk+1;
    subplot(3,3,kk)
    hold on
    yyaxis left

    snifft = interp1((0:length(sniff_trace(ss,:))-1)/10,sniff_trace(ss,:),timeaxis);
    plot(timeaxis,snifft)
    ylabel('Sniff Trace')

    yyaxis right
    % plot((1:wind)/10,mean(corrmod,1))
    pltter = squeeze(rsa_P1(ss,nn,:,2))-squeeze(rsa_P1(ss,nn,:,1));
    plot(timeaxis,pltter)
   

    ylabel('RSA performance')
    if and(ss==3,nn==3)
        legend({'Sniff','Mod-Unmod'})
    end
    xlabel('time(s)')
    title(sprintf('RSA for sub %02d, %s',ss,anat_names{nn}))    
end
end


rsa_wt = squeeze(sum(abs(rsa_P1),3));
ARC_barplot(rsa_wt)
xticks(1:nanat)
xticklabels(anat_names)
ylabel('Modulation avg')
% % mean plot
% figure('Position',[0 0 1480 720])
% kk = 0;
% for nn = 1:3
%     kk = kk+1;
%     subplot(3,1,kk)
%     hold on
%     yyaxis left
% 
%     snifft = interp1((0:length(sniff_trace(1,:))-1)/10,mean(sniff_trace),timeaxis);
%     plot(timeaxis,snifft)
%     ylabel('Sniff Trace')
% 
%     yyaxis right
%     % plot((1:wind)/10,mean(corrmod,1))
%     pltter = squeeze(rsa_P1(:,nn,:,2))-squeeze(rsa_P1(:,nn,:,1));
%     pltter = squeeze(mean(pltter,1));
%     plot(timeaxis,pltter)
% 
% 
%     ylabel('RSA performance')
%     if and(ss==3,nn==3)
%         legend({'Sniff','Mod-Unmod'})
%     end
%     xlabel('time(s)')
%     title(sprintf('RSA for sub %02d, %s',ss,anat_names{nn}))    
% end
% 

%% Figure
% rsa_pvals = tcdf(rsa_P1,sum(utl_mask,'all'),'upper');
% 
% load('C:\Work\SFP\common_datafiles\snifflabels.mat');
% num_desc = 8;
% narea = 3;
% figure()
% hold on
% rsa_P1(2,1,:,1) = nan;
% for nn = 1:narea
%     if nn==1
%         ylabel('Beta weight')
%     end
%     subplot(1,narea,nn)
%     hold on
%     ARC_barplot(abs(squeeze(rsa_P1(:,nn,1:num_desc,:))),false)
%     xticks(1:num_desc)
%     xticklabels(proper_list(1:num_desc))
%     xtickangle(90)
%     title(sprintf('Anat: %s',anat_names{nn}))
% end
% 
% S_mat = squeeze(mean(rsa_P1));
% rsa_pvals_mean = tcdf(S_mat,sum(utl_mask,'all'),'upper');
% 
% savefig(fullfile(savepath,'feat_map_comp'))
% print(fullfile(savepath,'feat_map_comp'),'-dpng')
% clearvars fless_mat mainmat Fless_mat_pruned Fless_mat unity task_run set_run sess_run anat_cell
% 
% % clear Fmat_1_m behav_mat unity M_anat M_reg n_reg M_unstack fless_mat fless_mat_unn modelmd_ modelmd S_omat_vals utl_mask utl_mask2
% save(fullfile(savepath,'SFP_snifffeat'),'settings_','rsa_P1')
