%% RSA - Trial-wise
numpcs = [14 11 11]; % 90% Variance
root = 'C:\Work';
corrmat_ = true;
if corrmat_
    settings_.nodor = 160;
    settings_.wind = 7500; % Number of samples
    settings_.nsniffcomp = 31;
    settings_.loadvec = [3 4 9:settings_.nsniffcomp];
    settings_.featorflessnot = true;
    settings_.chem = false;
    if  settings_.chem; delfeat = 1; else; delfeat = 0; end
    rsa_pvals = zeros(3,4);

    dirs = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
        fullfile(root ,'\SFP\sfp_behav_s02_correct');
        fullfile(root ,'\SFP\sfp_behav_s04_correct')};

    dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
        fullfile(root,'ARC\ARC\ARC02\single');
        fullfile(root,'ARC\ARC\ARC03\single')};

    behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));
    savepath = 'C:\Work\SFP\Final_plots\Behavioral\Trialwise RSA\temp';
    hold on
    nconds = 4;
    rsa_P1 = zeros(3,nconds,2+delfeat);
    rsa_P1t = zeros(3,nconds,2+delfeat);
    for ss = 1:length(dirs)
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
        Feat_mat_pruned(isnan(Feat_mat_pruned))=0;
        Feat_mat_pruned = zscore(Feat_mat_pruned,1);

        [coeff,Feat_mat_pruned,~,~,var] = pca(Feat_mat_pruned);
        Feat_mat_pruned = Feat_mat_pruned(:,1:numpcs(ss));


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
        utl_mask = logical(triu(ones(length(group_vec)),1)); % All possible odors

        % Behavioral RSMs
        behav_ratings = behav.behav(ss).ratings(group_vec,:);

        if settings_.featorflessnot; mainmat = Feat_mat_pruned; else; mainmat = Fless_mat_pruned; end
        mainmat(isnan(mainmat))=0;

        if settings_.featorflessnot
            mainmat = zscore(mainmat,[],1);
            % A2_corr = corrcoef(mainmat');
        else
            % A2_corr = pdist(mainmat,"spearman");
            % A2_corr = 1-squareform(A2_corr);
        end
        A2_corr = pdist(mainmat,"correlation");
        A2_corr = 1-squareform(A2_corr);

        behav_corr = corrcoef(behav_ratings');
        int_corr = -abs(behav_ratings(:,1)-behav_ratings(:,1)');
        pls_corr = -abs(behav_ratings(:,2)-behav_ratings(:,2)');


        [wt,t_sc] = ARC_multicomputeWeights_tsc([behav_corr(utl_mask) unity(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)], A2_corr(utl_mask) );
        if settings_.chem
            load(fullfile(statpath,'chem_corr.mat'))
            chem_corr = Behav_RSM(group_vec,group_vec);
            [wt,t_sc] = ARC_multicomputeWeights_tsc([behav_corr(utl_mask) unity(utl_mask) chem_corr(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)], A2_corr(utl_mask) );
        end

        
        
        rsa_P1(ss,1,1:2+delfeat) = wt(2:3+delfeat);
        rsa_P1t(ss,1,1:2+delfeat) = t_sc(2:3+delfeat);

        
        rsa_pvals(ss,1)=  SFP_crosspvalcalc(wt(2:3+delfeat),t_sc(2:3+delfeat),sum(utl_mask(:)));

        if settings_.chem
            [wt,t_sc] = ARC_multicomputeWeights_tsc([behav_corr(utl_mask) unity(utl_mask) chem_corr(utl_mask) int_corr(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)],  A2_corr(utl_mask));
        else
            [wt,t_sc] = ARC_multicomputeWeights_tsc([behav_corr(utl_mask) unity(utl_mask) int_corr(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)],  A2_corr(utl_mask));

        end
        
        rsa_P1(ss,2,1:2+delfeat) = wt(2:3+delfeat);
        rsa_P1t(ss,2,1:2+delfeat) = t_sc(2:3+delfeat);

        rsa_pvals(ss,2)=  SFP_crosspvalcalc(wt(2:3+delfeat),t_sc(2:3+delfeat),sum(utl_mask(:)));

        
        if settings_.chem
            [wt,t_sc] = ARC_multicomputeWeights_tsc([behav_corr(utl_mask) unity(utl_mask) chem_corr(utl_mask) pls_corr(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)],  A2_corr(utl_mask));
        else
            [wt,t_sc] = ARC_multicomputeWeights_tsc([behav_corr(utl_mask) unity(utl_mask) pls_corr(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)],  A2_corr(utl_mask));
        end

        rsa_P1(ss,3,1:2+delfeat) = wt(2:3+delfeat);
        rsa_P1t(ss,3,1:2+delfeat) = t_sc(2:3+delfeat);

          rsa_pvals(ss,3)=  SFP_crosspvalcalc(wt(2:3+delfeat),t_sc(2:3+delfeat),sum(utl_mask(:)));

            
        if settings_.chem
            [wt,t_sc] = ARC_multicomputeWeights_tsc([behav_corr(utl_mask) unity(utl_mask) chem_corr(utl_mask)  int_corr(utl_mask)  pls_corr(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)],  A2_corr(utl_mask));
        else
            [wt,t_sc] = ARC_multicomputeWeights_tsc([behav_corr(utl_mask) unity(utl_mask)  int_corr(utl_mask)  pls_corr(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)],  A2_corr(utl_mask));
        end
        rsa_P1(ss,4,1:2+delfeat) = wt(2:3+delfeat);
        rsa_P1t(ss,4,1:2+delfeat) = t_sc(2:3+delfeat);

          rsa_pvals(ss,4)=  SFP_crosspvalcalc(wt(2:3+delfeat),t_sc(2:3+delfeat),sum(utl_mask(:)));
    end
    p_values_3d = ARC_RSA_pvals(rsa_P1t, rsa_P1, sum(utl_mask(:)))
     % p_values_3d = ARC_RSA_pvals(rsa_P1t(1,:,:), rsa_P1(1,:,:), sum(utl_mask(:)))


    ARC_barplot(rsa_P1)
    % rsa_P1 = rsa_P1t;
    xticks(1:4)
    xticklabels({'Sniff RSA','-Int','-Pls','-Int-Pls'})
    ylabel('Representational Similarity (std. \beta)')
    legend({'Perceptual similarity','Odor trial similarity'})

    % yline(r2t(0.05,sum(utl_mask2(:))));
    % yline(r2t(0.05,nchoosek(length( group_vec),2)));
    savefig(fullfile(savepath,'fless_map'))
    print(fullfile(savepath,'fless_map'),'-dpng')
    SFP_clearLargeVariables(10)
    % clear Fmat_1_m behav_mat Fless_mat Fless_mat_pruned A2_corr behav_corr int_corr pls_corr sess_run set_run task_run unity M_anat M_reg n_reg M_unstack fless_mat fless_mat_unn modelmd_ modelmd S_omat_vals utl_mask utl_mask2
    save(fullfile(savepath,'ARC_RSA_fless'))
end