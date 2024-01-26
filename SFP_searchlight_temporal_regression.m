%% General Settings
root = 'C:\Work';
nodor = 160;
twind = 75; % Number of samples

tidx = 1:5:65;
wind = length(tidx);
idx_delay = [5 7 5 5];

dirs = {fullfile(root ,'\SFP\sfp_behav_s01_peak');
        fullfile(root ,'\SFP\sfp_behav_s02_peak');
        fullfile(root ,'\SFP\sfp_behav_s03_peak')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
        fullfile(root,'ARC\ARC\ARC02\single');
        fullfile(root,'ARC\ARC\ARC03\single')};

dirs3 = {fullfile(root,'SFP\Flow_temporal_regress\SFP1');
        fullfile(root,'SFP\Flow_temporal_regress\SFP2');
        fullfile(root,'SFP\Flow_temporal_regress\SFP3')};

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
rsa_P1_ = cell(3,nanat); % Complete
rsa_P2_ = cell(3,nanat); % Complete

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

    load(fullfile(statpath,'sfp_feats_corr.mat'))
    Fless_mat = vertcat(fless_mat_peak{:});
    Fless_mat_pruned = Fless_mat(:,1:twind);

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
    tnvox = sum(anatmasks,'all');
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
    h = waitbar(0, 'Please wait...'); % Progress bar per subject
    k = 0;
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
        
        rsa_vec = zeros(nvox,wind);
        rsa_vec2 = zeros(nvox,wind);
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
                for tpoint = tidx
                    k = k+1;
                    waitbar(k/(tnvox*wind), h, sprintf('Processing S:%02d, %.2f %%', ss, (k/(tnvox*wind)*100))) % Update progress

                    Fmat = Fless_mat_pruned(:,tpoint);
                    Fmat_m = -abs(Fmat-Fmat');

                    Fmat2 = Fless_mat_pruned(:,tpoint+idx_delay(ss));
                    Fmat_m2 = -abs(Fmat2-Fmat2');
    
                    temp = SFP_computeWeights(Fmat_m(utl_mask),Fmat_m2(utl_mask),M_reg);
                    rsa_vec(cnt2,tpoint) = temp(2);
                    rsa_vec2(cnt2,tpoint) = temp(3);
                end
            end
        end
        rsa_P1_{ss,ii} = rsa_vec;
        rsa_P2_{ss,ii} = rsa_vec2;

        % if mapper
        %     for tpoint = 1:wind
        %         rsa_vec_3d = unmasker(rsa_vec(:,tpoint),logical(anatmasks(:,:,:,ii)));
        %         write_reshaped_nifty(rsa_vec_3d, savepath, false, fullfile(anatpath,maskfile), sprintf('SFP%02d_%s_t%02d',ss,anat_names{ii},tpoint));
        %     end
        % end
    end
    close(h);
end
clear Fmat_1_m Fmat_2_m unity M_anat M_reg n_reg M_unstack fless_mat fless_mat_unn modelmd_ modelmd S_omat_vals utl_mask utl_mask2
save(fullfile(savepath,'ARC_RSA'))
% 
% Group Plot
load('C:\Work\SFP\RSA_timecourse\timecourse_peak.mat')
neural_RSA = cellfun(@mean,rsa_P1_,'UniformOutput',false); 

neural_RSA2 = cellfun(@mean,rsa_P2_,'UniformOutput',false); 

figure('Position',[0 0 1480 220])
kk = 0;
for nn = 1:nanat
    kk = kk+1;
    subplot(1,nanat,kk)
    hold on
    yyaxis left
    plot((1:twind)/10,mean(sniff_trace,1))
    ylabel('Sniff Trace')

    yyaxis right
    % plot((1:wind)/10,mean(corrmod,1))
    t1 = mean(vertcat(neural_RSA{:,nn}));
    plot((tidx)/10, t1(tidx))
    t2 = mean(vertcat(neural_RSA2{:,nn}));
    plot((tidx)/10,  t2(tidx))

    ylabel('RSA performance')
    % legend({'Sniff','Neural Similarity'})
    xlabel('time(s)')
    if nn==nanat
    legend({'Sniff','Sniff-RSA','Perceptual-RSA'})
    end
    title(sprintf('RSA for %s',anat_names{nn}))    
end
savefig(fullfile(savepath,'ARC_RSA'))
print(fullfile(savepath,'ARC_RSA'),'-dpng')

% Subject-wise plot
figure('Position',[0 0 1480 720])
kk = 0;
for ss = 1:3
for nn = 1:nanat
    kk = kk+1;
    subplot(3,nanat,kk)
    hold on
    yyaxis left
    plot((1:twind)/10,sniff_trace(ss,:))
    ylabel('Sniff Trace')

    yyaxis right
    % plot((1:wind)/10,mean(corrmod,1))
    t1 = neural_RSA{ss,nn};
    plot((tidx)/10, t1(tidx))
    t2 = neural_RSA2{ss,nn};
    plot((tidx)/10,  t2(tidx))


    ylabel('RSA performance')
    % legend({'Sniff','Neural Similarity'})
    xlabel('time(s)')
    title(sprintf('RSA for sub %02d, %s',ss,anat_names{nn}))    
end
end
savefig(fullfile(savepath,'ARC_RSA_sub'))
print(fullfile(savepath,'ARC_RSA_sub'),'-dpng')

% % rsa_P1_pvals = cellfun(@(x) r2p(x,nchoosek(ntrials,2)),rsa_P1_,'UniformOutput',false);
% % rsa_P1_pvals2 = cat(1,rsa_P1_pvals{:});
% % rsa_P1_thresh = fdr_benjhoc(rsa_P1_pvals2);
% 
% % for pp = 1:numel(rsa_P1_thresh);  if isempty(rsa_P1_thresh{pp}); rsa_P1_thresh{pp} = 0; end; end
% % rsa_P1_thresh{1,1,1} = 0; % Check this manually
% % rsa_P1_thresh{3,1,2} = 0; % Check this manually
% 
% % rsa_P1 = cellfun(@(x) (sum(x>r2t(0.05,nchoosek(ntrials,2)))./length(x))*100,rsa_P1_);
% % rsa_P1 = cellfun(@(x1,x2) (sum(x1<x2.*ones(size(x1)))./length(x1))*100,rsa_P1_pvals,rsa_P1_thresh);
% % rsa_P1 = cellfun(@(x) (sum(x<rsa_P1_thresh)./length(x))*100,rsa_P1_pvals);
% rsa_P1 = cellfun(@(x) mean(x),rsa_P1_);
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
% legend('peak','rsa')
% ylabel('Representational Similarity (r)')
% % yline(r2t(0.05,sum(utl_mask2(:))));
% yline(r2t(0.05,nchoosek(length( group_vec),2)));
% savefig(fullfile(savepath,'ARC_RSA'))
% print(fullfile(savepath,'ARC_RSA'),'-dpng')
