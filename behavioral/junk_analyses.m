%% Correlation with behavior
% Sniffing modulation
grp_ = true;
rootf = 'C:\Work\SFP\Final_plots\Behavioral';
if grp_
    nodor = 160;
    wind = 3500; % Number of samples
    dirs = {'C:\Work\SFP\sfp_behav_s01_correct';
        'C:\Work\SFP\sfp_behav_s02_correct';
        'C:\Work\SFP\sfp_behav_s04_correct'};
    color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];

    behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));
    behav_ids = [1 1 1; 2 2 2; 3 3 3; 17 15 15; 8 8 8; 9 9 9];
    behav_names = behav.behav(2).percepts(behav_ids(:,2));


    figure('Position',[0 0 1280 720])
    hold on
    plot_id = 0;
    for ss = 1:length(dirs)
        load(fullfile(dirs{ss},'sfp_feats_corrected.mat'))
        Fless_mat = vertcat(fless_mat{:});
        feat_mat = vertcat(feat_mat{:});
        feat_mat_pruned = feat_mat(:,[3 4 9:14]);

        anatdir = fullfile('C:\Work\ARC\ARC\',sprintf('ARC%02d',ss),'single');

        if ss==3; s2 = 4; else; s2 = ss; end
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

        for pp = 1:size(behav_ids,1)
            plot_id = plot_id +1;
            v_id = behav_ids(pp,ss);
            if v_id>0
                behav_ratings = behav.behav(ss).ratings(:,v_id);
                behav_ratings_= behav_ratings(group_vec);
                divs = quantile(behav_ratings_,[0.1 0.5 0.9]);
                divs2 = quantile(behav_ratings_,19);
                [~,argsort] = sort(behav_ratings_);

                subplot(3,size(behav_ids,1),plot_id)
                hold on
                xt = 0:0.001:(wind-1)/1000;
                % v_neg_hard = Fless_mat(behav_ratings_<= divs(1),1:wind );
                % v_pos_hard = Fless_mat(behav_ratings_>= divs(3),1:wind );

                v_neg_hard = Fless_mat(behav_ratings_<= divs2(2),1:wind );
                v_pos_hard = Fless_mat(behav_ratings_>= divs2(end-1),1:wind );
                v_null_hard = Fless_mat(and(behav_ratings_>= divs2(9),behav_ratings_<= divs2(11)),1:wind );


                feat_neg_hard = feat_mat(behav_ratings_<= divs(1),[3 4 9:14]);
                feat_pos_hard = feat_mat(behav_ratings_>= divs(3),[3 4 9:14]);
                % 
                shaded_plot(xt,mean(v_neg_hard,1),1.96*std(v_neg_hard,1)/sqrt(size(v_neg_hard,1)),color(1,:));
                % %             shaded_plot(xt,mean(v_neg_soft,1),1.96*std(v_neg_soft,1)/sqrt(size(v_neg_soft,1)),color(2,:));
                % %             shaded_plot(xt,mean(v_pos_soft,1),1.96*std(v_pos_soft,1)/sqrt(size(v_pos_soft,1)),color(3,:));
                % shaded_plot(xt,mean(v_null_hard,1),1.96*std(v_null_hard,1)/sqrt(size(v_null_hard,1)),color(2,:));
                shaded_plot(xt,mean(v_pos_hard,1),1.96*std(v_pos_hard,1)/sqrt(size(v_pos_hard,1)),color(4,:));

                % imagesc(Fless_mat(argsort,1:wind))
                axis tight
                colormap('jet')
                % 
                % tic
                % y = [zeros(size(feat_neg_hard,1),1); ones(size(feat_pos_hard,1),1)];
                % g = cat(1,feat_neg_hard,feat_pos_hard);
                % g = mat2cell(g,size(g,1),ones(size(g,2),1));
                % p =  anovan(y,g);
                % toc
                if plot_id==3*size(behav_ids,1)
                    legend({'low','high'})
                end
                if ss ==1 
                title(sprintf('%s',behav_names{pp}))
                end
            end
            if and(ss==3,pp==3); xlabel('Time(s)'); end
            if and(ss==2,pp==1); ylabel('Amplitude'); end
        end
    end
    savefig(fullfile(rootf,'snifftraces_mat'))
    print(fullfile(rootf,'snifftraces_mat'),'-dpng')
end

%% Behavioral quantification
% Sniffing modulation
tic
grp_ = true;
rootf = 'C:\Work\SFP\Final_plots\Behavioral';
if grp_
    nodor = 160;
    wind = 3500; % Number of samples
    dirs = {'C:\Work\SFP\sfp_behav_s01_correct';
        'C:\Work\SFP\sfp_behav_s02_correct';
        'C:\Work\SFP\sfp_behav_s04_correct'};
    color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
    behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));
    ndisc = size(behav.behav(1).ratings,2);

    num_eff = zeros(3,ndisc);
    num_pval = zeros(3,ndisc);
    for ss = 1:length(dirs)
        fprintf('subject %02d: \n',ss)
        load(fullfile(dirs{ss},'sfp_feats_corrected.mat'))
        Fless_mat = vertcat(fless_mat{:});
        feat_mat = vertcat(feat_mat{:});
        feat_mat_pruned = feat_mat(:,[3 4 9:14]);

        anatdir = fullfile('C:\Work\ARC\ARC\',sprintf('ARC%02d',ss),'single');

        if ss==3; s2 = 4; else; s2 = ss; end
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

        for pp = 1: ndisc

      
            behav_ratings = behav.behav(ss).ratings(:,pp);
            behav_ratings_= behav_ratings(group_vec);
           
            nbins = 3;
            delta = 0.0001*linspace(0,1,nbins+1);
            edges_ = quantile(behav_ratings_, linspace(0, 1, nbins + 1));
            edges_ = edges_+delta;
            y = discretize(behav_ratings_, edges_);
            
            [ num_eff(ss,pp) , num_pval(ss,pp) ] = SFP_logisticRegressionCV(y, feat_mat_pruned);
    
        end

    end
end 

nS = 3;
behav_lab1 = behav.behav(1).percepts;
behav_lab2 = behav.behav(2).percepts;

% for ii=1:nS
%     behav(ii).rel(end+1)=0;
% end

idx2 = [1:10 19 11:14 19 15 19 16:18]; % Use this is argsort in sub(2 and 3) to match the labels
idx1 = [1:18 19 19 19];
behav_labs = {behav_lab1{:} behav_lab2{end-2:end}};

num_eff(:,end+1)=nan;
idx1 = [1:18 19 19 19];

rels = zeros(nS,21);
for ii = 1:3
    if ii==1
        rels(ii,:)=num_eff(ii,idx1);
    else       
        rels(ii,:)=num_eff(ii,idx2);
    end
end
rels_mean = mean(rels,1);
[rels_sort,argsort] = sort(rels_mean,'descend');
behav_labs_sort = behav_labs(argsort);

figure()
bar(1:21,nanmean(rels))
hold on
errorbar(1:21,nanmean(rels),1.96*nanstd(rels)/sqrt(3),'.')
c = {'r.','g.','b.'};
for ss = 1:3
    plot(1:21,rels(ss,:),c{ss})
end
xticks(1:21)
xticklabels(behav_labs)
xtickangle(90)
yline(1/3)

toc

%% RSA
root = 'C:\Work';
corrmat_ = true;
numpcs = [14 11 11]; % 90% Variance
if corrmat_
    settings_.nodor = 160;
    settings_.wind = 7500; % Number of samples
    settings_.nsniffcomp = 32;
    settings_.loadvec = [3 4 9:settings_.nsniffcomp];
    settings_.featorflessnot = true;

    dirs = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
        fullfile(root ,'\SFP\sfp_behav_s02_correct');
        fullfile(root ,'\SFP\sfp_behav_s04_correct')};

    dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
        fullfile(root,'ARC\ARC\ARC02\single');
        fullfile(root,'ARC\ARC\ARC03\single')};

    behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));
    savepath = 'C:\Work\SFP\Final_plots\Basic RSA';
    hold on
    nconds = 4;
    rsa_P1 = zeros(3,nconds);
    for ss = 1:length(dirs)
        fprintf('Subject: %02d\n',ss)
        if ss==3; s2 = 4; else; s2 = ss; end
        statpath = dirs{ss};
        anatdir = dirs2{ss};
        % savepath = dirs3{ss};
        mkdir(savepath)

        load(fullfile(statpath,'sfp_feats_main.mat'))
        load(fullfile(statpath,'task_struct.mat'))
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
        utl_mask = logical(triu(ones(length(unique(group_vec))),1)); % All possible odors

        % Behavioral RSMs
        behav_ratings = behav.behav(ss).ratings;

        if settings_.featorflessnot; mainmat = Feat_mat_pruned; else; mainmat = Fless_mat_pruned; end
        mainmat(isnan(mainmat))=0;

        if settings_.featorflessnot
            mainmat = zscore(mainmat,[],1);
        end
        
        A2 = splitapply(@mean,mainmat,group_vec);
        A2_corr = corrcoef(A2');

        behav_corr = corrcoef(behav_ratings');
        int_corr = -abs(behav_ratings(:,1)-behav_ratings(:,1)');
        pls_corr = -abs(behav_ratings(:,2)-behav_ratings(:,2)');


        [~,t_sc] = ARC_multicomputeWeights_tsc([A2_corr(utl_mask)  task_struct(utl_mask) task_run(utl_mask)], behav_corr(utl_mask));
        rsa_P1(ss,1) = t_sc(2);
     
        [~,t_sc] = ARC_multicomputeWeights_tsc([A2_corr(utl_mask) int_corr(utl_mask) task_struct(utl_mask) task_run(utl_mask)], behav_corr(utl_mask));
        rsa_P1(ss,2) = t_sc(2);

       [~,t_sc] = ARC_multicomputeWeights_tsc([A2_corr(utl_mask) pls_corr(utl_mask) task_struct(utl_mask) task_run(utl_mask)], behav_corr(utl_mask));
        rsa_P1(ss,3) = t_sc(2);

        [~,t_sc] = ARC_multicomputeWeights_tsc([A2_corr(utl_mask) int_corr(utl_mask)  pls_corr(utl_mask) task_struct(utl_mask) task_run(utl_mask)], behav_corr(utl_mask));
        rsa_P1(ss,4) = t_sc(2);

    end
    S_mat = squeeze(mean(rsa_P1));
    S_err = squeeze(std(rsa_P1))./sqrt(3);
    figure('Position',[0.5 0.5 400 250])
    hold on
    ngroups = size(S_mat, 1);
    nbars = size(S_mat, 2);
    bar(1:nconds,S_mat);
    errorbar(1:nconds,S_mat,S_err,'.')
    yline(tinv(0.95,sum(utl_mask(:))))
 
    % legend({'Perceptual','Chemical','Mutual'})
    % legend()
    % Subject data points
    c_s = {'r','g','b'}; % Data dots for subjects
        for jj = 1:3
            plot([1:nconds],squeeze(rsa_P1(jj,:)),c_s{jj},'handle','off')
        end
    % legend('S_I','S_O','sniff col.','neural collinearity')
    xticks(1:4)
    xticklabels({'Sniff RSA','-Int','-Pls','-Int-Pls'})
    ylabel('Representational Similarity (t)')
    % yline(r2t(0.05,sum(utl_mask2(:))));
    % yline(r2t(0.05,nchoosek(length( group_vec),2)));
    savefig(fullfile(savepath,'feat_map'))
    print(fullfile(savepath,'feat_map'),'-dpng')

    % clear Fmat_1_m behav_mat unity M_anat M_reg n_reg M_unstack fless_mat fless_mat_unn modelmd_ modelmd S_omat_vals utl_mask utl_mask2
    save(fullfile(savepath,'ARC_RSA'),'settings_','rsa_P1')
end

%% Raw Sniff plots

% Use ss = 2 and extract group_vec, Fless_mat_pruned for s 2
s1s = Fless_mat_pruned(group_vec==157,:);
s2s = Fless_mat_pruned(group_vec==142,:);
% plot(s1s','r')
% hold on
% plot(s2s','b')

figure()
shaded_plot(1:size(s1s,2),mean(s1s),std(s1s)/sqrt(size(s1s,1)),'r')
hold on
shaded_plot(1:size(s1s,2),mean(s2s),std(s2s)/sqrt(size(s2s,1)),'b')
xticks([0 2500 5000 7500])
xticklabels([0 2.5 5 7.5])
xlabel('Time(s)')
ylabel('Sniff flow (AU)')
legend('3-(5-Methyl-2-furyl)butanal','Ethyl 2-methylpentanoate')


%% RSA time-series
corrmat_ = true;
dwnsample = 100;
if corrmat_
    nodor = 160;
    wind = 75; % Number of samples
    dirs = {'C:\Work\SFP\sfp_behav_s01_correct';
        'C:\Work\SFP\sfp_behav_s02_correct';
        'C:\Work\SFP\sfp_behav_s04_correct'};
    %   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
    behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));

    hold on
    corrmod = zeros(3,wind);
    corrmod_int = zeros(3,wind);
    corrmod_pls = zeros(3,wind);
    sniff_trace = zeros(3,wind);
    for ss = 1:length(dirs)
        load(fullfile(dirs{ss},'sfp_feats_main.mat'))
        Fless_mat = vertcat(fless_mat{:});
        anatdir = fullfile('C:\Work\ARC\ARC\',sprintf('ARC%02d',ss),'single');

        if ss==3; s2 = 4; else; s2 = ss; end
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

        Fless_mat_pruned = Fless_mat(:,1:dwnsample:end);

        % Behavioral features
        behav_ratings = behav.behav(ss).ratings;
        behav_ratings_ =  behav_ratings(group_vec,:);
        b1 = corrcoef(behav_ratings_')';

        rint = regressmeout(behav_ratings(:,[2:end])',repmat(behav_ratings(:,1)',17,1));
        rint = rint(:,group_vec);
        brint = corrcoef(rint);

        rpls = regressmeout(behav_ratings(:,[1 3:end])',repmat(behav_ratings(:,2)',17,1));
        rpls = rpls(:,group_vec);
        bpls = corrcoef(rpls);

        for jj = 1:wind
            Fless_corr = -abs(Fless_mat_pruned(:,jj)-Fless_mat_pruned(:,jj)');
            b1_vec = b1(utl_mask);
            b1_vec = regressmeout(b1_vec',unity(utl_mask)')';

            temp = corrcoef(b1_vec ,Fless_corr( utl_mask ));
            corrmod(ss,jj) = temp(2);

            % temp = corrcoef(brint,Fless_corr);
            % corrmod_int(ss,jj) = temp(2);

            % temp = corrcoef(bpls,Fless_corr);
            % corrmod_pls(ss,jj) = temp(2);
        end
        sniff_trace(ss,:) = mean(Fless_mat_pruned ,1);
    end
end

figure('Position',[0.5 0.5 640 480])
hold on
yyaxis left
plot((1:wind)/10,mean(sniff_trace,1))
ylabel('Sniff Trace')

yyaxis right
plot((1:wind)/10,mean(corrmod,1))
% plot((1:wind)/10,mean(corrmod_int,1))
% plot((1:wind)/10,mean(corrmod_pls,1))
ylabel('RSA performance')

legend({'Sniff','Perceptual Similarity'})
yyaxis left
xlabel('time(s)')

figure()
hold on
for ii = 1:3
    subplot(1,3,ii)
    hold on
    yyaxis left
    plot((1:wind)/10,corrmod(ii,:))
    plot((1:wind)/10,corrmod_int(ii,:))
    plot((1:wind)/10,corrmod_pls(ii,:))
    ylabel('RSA performance')
    yyaxis right
    plot((1:wind)/10,sniff_trace(ii,:))
    ylabel('Sniff Trace')
    legend({'Basic','-Int','-Pls','Sniff'})
    yyaxis left
    title(sprintf('Subject: %02d',ii))
end

figure()
hold on
x = mean(sniff_trace,1);
y = mean(corrmod,1);
c = linspace(1,10,wind);
scatter(x(1:75),y(1:75),25,c(1:75))
ylabel('Rep Similarity')
xlabel('Sniff')
% 
figure('Position',[0 0 720 200])
hold on
for ii = 1:3
    subplot(1,3,ii)
    hold on
    x = sniff_trace(ii,:);
    y = corrmod(ii,:);
    c = linspace(1,10,wind);
    scatter(x,y,25,c)
    title(sprintf('Subject: %02d',ii))
end
savefig('pls_sub_cross')
print('pls_sub_cross','-dpng')

%% RSA static bootstrap
root = 'C:\Work';

ncomps = [2 3 4 4];
nodor = 160;
wind = 75; % Number of samples
dirs = {fullfile(root ,'\SFP\sfp_behav_s01');
        fullfile(root ,'\SFP\sfp_behav_s02');
        fullfile(root ,'\SFP\sfp_behav_s03')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
        fullfile(root,'ARC\ARC\ARC02\single');
        fullfile(root,'ARC\ARC\ARC03\single')};

dirs3 = {fullfile(root,'SFP\CCA_2comp\SFP1');
        fullfile(root,'SFP\CCA_2comp\SFP2');
        fullfile(root,'SFP\CCA_2comp\SFP3')};

%   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));
hold on
nboot = 10;
corrmod = zeros(3,nboot);

for ss = 1:length(dirs)
    statpath = dirs{ss};
    anatdir = dirs2{ss};
    savepath = dirs3{ss};


    onsets = load(fullfile(anatdir,sprintf('conditions_NEMO%02d.mat',ss)),'onsets');
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

    behav_ratings = behav.behav(ss).ratings;
    behav_ratings = behav_ratings(group_vec,:);


    % Sniff descriptors
    load(fullfile(dirs{ss},'sfp_feats.mat'))
    fless_mat{1,19}(:,701:end) = [];
    Fless_mat = vertcat(fless_mat{:});
    % Fless_mat_pruned = Fless_mat(:,1:wind);
    
    Feat_mat_pruned = vertcat(feat_mat{:});
    Feat_mat_pruned =  Feat_mat_pruned(:,[3 4 9:14]) ;

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
    Fless_mat_pruned = U;

    anatdir = fullfile('C:\Work\ARC\ARC\',sprintf('ARC%02d',ss),'single');

    if ss==3; s2 = 4; else; s2 = ss; end
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


    Fless_corr = corrcoef(Feat_mat_pruned');
    Fless_corr(isnan(Fless_corr))=0;

    % Behavioral features
    behav_ratings = behav.behav(ss).ratings;
    behav_ratings_ =  behav_ratings(group_vec,:);

    b1 = corrcoef(behav_ratings_')';
    b1_vec = b1(utl_mask);
    b1_vec = regressmeout(b1_vec',unity(utl_mask)')';
    b1 = SFP_vector_to_matrix(b1_vec,size(b1));
    [corrmod(ss,:),p] = SFP_bootstrap_corr(b1,Fless_corr,nboot);
end

figure()
 m = mean(corrmod,2);
bar(mean(m))
hold on
errorbar([1],mean(m),std(m),'.')
plot([1,1,1],m,'k.','Markersize',10)
ylabel('representational similarity')
yline(r2t(0.05,nchoosek(4560,2)))
xticks([1])
xticklabels({'perceptual RSA with sniff'})



%% RSA static - CCA
cca_maker = true;
if cca_maker
root = 'C:\Work';
autoreducer = false;
presniff_clear = true;
ncomps = [5 5 5 5];
% ncomps = [2 3 4 4 ];
nodor = 160;
wind = 75; % Number of samples
dirs = {fullfile(root ,'\SFP\sfp_behav_s01_fullfeat_14');
        fullfile(root ,'\SFP\sfp_behav_s02_fullfeat_14');
        fullfile(root ,'\SFP\sfp_behav_s04_fullfeat_14')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
        fullfile(root,'ARC\ARC\ARC02\single');
        fullfile(root,'ARC\ARC\ARC03\single')};

dirs3 = {fullfile(root,'SFP\CCA_2comp_full\SFP1');
        fullfile(root,'SFP\CCA_2comp_full\SFP2');
        fullfile(root,'SFP\CCA_2comp_full\SFP3')};

nplt = 10;
%   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));
hold on

corrmod = zeros(3,2);
Amat = {};
figure1 = figure();
figure2 = figure();

  sbplt  =0;
for ss = 1:length(dirs)
  
    statpath = dirs{ss};
    anatdir = dirs2{ss};
    savepath = dirs3{ss};

    % Perceptual descriptors
    onsets = load(fullfile(anatdir,sprintf('conditions_NEMO%02d.mat',ss)),'onsets');
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
    utl_mask = logical(triu(ones(length(unity)),1)); % All possible odors

    behav_ratings = behav.behav(ss).ratings;
    behav_ratings = behav_ratings(group_vec,:);

    % Sniff descriptors
    load(fullfile(dirs{ss},'sfp_feats_main.mat'))

  

    fless_mat{1,19}(:,701:end) = [];
    % Fless_mat = vertcat(fless_mat{:});
    % Fless_mat_pruned = Fless_mat(:,1:wind);
    
    Feat_mat_pruned = vertcat(feat_mat{:});
    Feat_mat_pruned =  Feat_mat_pruned(:,[3 4 9:end]) ;

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


    if presniff_clear
        behav_ratings = behav_ratings(sniff_cmat,:);
        utl_mask = utl_mask(sniff_cmat,sniff_cmat);
        unity = unity(sniff_cmat,sniff_cmat);
        F1_no_near_constant = F1_no_near_constant(sniff_cmat,:);
        F2_no_near_constant = F2_no_near_constant(sniff_cmat,:);
    end

    [F1_coeff, F2_coeff, r, U,V,stats] = canoncorr(F1_no_near_constant, F2_no_near_constant);

    figure(figure1);
     for plt = 1:nplt
        sbplt = sbplt+1;
        subplot(3,nplt,sbplt) 
        bar(F1_coeff(:,plt))    
        xticklabels({'Inflow','ExFlow','Time Peak','Time Trough','InVol','ExVol','InDur','ExDur'})
        title(sprintf('CCA%01d r:%.2f',plt,r(plt)))
    end

    figure(figure2);
    subplot(1,3,ss)
    plot(r)

    % Correlations
    Fless_corr = corrcoef(Feat_mat_pruned');
    Fless_corr(isnan(Fless_corr))=0;

    if autoreducer; ncomps_= sum(r>r2t(0.05,size(F1_no_near_constant,1))); else; ncomps_ = ncomps(ss); end
    
    U = U(:,1:ncomps_);
    A1 = corrcoef( F1_no_near_constant');
    A2 = corrcoef(U');

    Amat{ss} = A2;
    % Behavioral features
    b1 = corrcoef(behav_ratings')';
    b1_vec = b1(utl_mask);
    b1_vec = regressmeout(b1_vec',unity(utl_mask)')';
    tmp = corrcoef(b1_vec,A1(utl_mask));
    corrmod(ss,1) = tmp(2);
    tmp  = corrcoef(b1_vec,A2(utl_mask));
    corrmod(ss,2) = tmp(2);  
end

figure()
m = mean(corrmod,1);
bar(m)
hold on
errorbar([1 2],m,std(corrmod),'.')
for ss = 1:3; plot([1,2],corrmod(ss,:),'k','Marker','.','Markersize',10); end
ylabel('Pattern decoding')
yline(r2t(0.05,nchoosek(4560,2)))
xticks([1 2])
xticklabels({'Basic','CCA'})
end


%% Make SVM - libsvm - Fmat  and behav
nfolds = 10;

if svm_trainer3
corrmod = zeros(nfolds,3);
corrmod2 = cell(3,1);
for ss = 1:length(dirs)
    load(fullfile(dirs{ss},'sfp_feats.mat'))
    behav_ratings = behav.behav(ss).ratings;
    Fless_mat = vertcat(feat_mat{:});
    anatdir = fullfile('C:\Data\ARC\',sprintf('ARC%02d',ss),'single');
    
    if ss==3; s2 = 4; else; s2 = ss; end
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
    
    Fless_mat_pruned = Fless_mat(:,[3 4 9:14]);
    Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
    oid_ = 1:160;
    oid = oid_(group_vec)';
    
    cvind = crossvalind('Kfold',size(Fless_mat_pruned,1),nfolds);
    pred_mat = zeros(size(Fless_mat_pruned,1)/nfolds,nfolds);
    ind_mat = zeros(size(Fless_mat_pruned,1)/nfolds,nfolds);
    for folderr = 1:nfolds
        ind = cvind==folderr;
        ind_mat(:,folderr) = find(ind); 
        X_train = Fless_mat_pruned(~ind,:);
        y_train = oid(~ind);
        X_test = Fless_mat_pruned(ind,:);
        y_test = oid(ind);
        mdl = svmtrain(y_train, X_train, ['-s 1 -t 2 -q']);
        [a,p] = svmpredict( y_test,  X_test, mdl,['-q']);

        corrmod(folderr,ss) = p(1);
        pred_rat = behav_ratings(a,:);
        test_rat = behav_ratings(y_test,:);
        pred_mat(:,folderr) = iter_corr(pred_rat,test_rat);      
    end
    corrmod2{ss} = pred_mat(:);
end
end

figure('Position',[0 0 320 240])
hold on
bar(mean(corrmod,1))
yline(100/160)
xticks(1:3)
xticklabels({'S1','S2','S3'})
ylabel('Performance (%)')
savefig('svm')
print('svm','-dpng')

figure('Position',[0 0 320 240])
corrmod2{1}(4321:end) = []; % This is wrong, but I am lazy
corr_m = horzcat(corrmod2{:});
boxplot(corr_m)
xticks(1:3)
xticklabels({'S1','S2','S3'})
ylabel('Performance (r)')
savefig('behav')
print('behav','-dpng')


p = [];
for ii = 1:3
     p(ii) = signrank(corrmod2{ss});
end
clear unity fless_mat fless_mat_unn Fless_mat utl_mask
save('svmpred')
%% Make SVM - libsvm - Fmat  - probabilistic
nfolds = 10;

if svm_trainer4
corrmod = zeros(nfolds,3);
corrmod2 = zeros(nfolds,3);
for ss = 1:length(dirs)
    load(fullfile(dirs{ss},'sfp_feats.mat'))
    Fless_mat = vertcat(feat_mat{:});
    anatdir = fullfile('C:\Data\ARC\',sprintf('ARC%02d',ss),'single');
    
    if ss==3; s2 = 4; else; s2 = ss; end
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
    
      Fless_mat_pruned = Fless_mat(:,[3 4 9:14]);
        Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
    oid_ = 1:160;
    oid = oid_(group_vec)';
    
    cvind = crossvalind('Kfold',size(Fless_mat_pruned,1),nfolds);
    pred_mat = zeros(size(Fless_mat_pruned,1)/nfolds,nfolds);
    ind_mat = zeros(size(Fless_mat_pruned,1)/nfolds,nfolds);
    for folderr = 1:nfolds
        ind = cvind==folderr;
        ind_mat(:,folderr) = find(ind); 
        X_train = Fless_mat_pruned(~ind,:);
        y_train = oid(~ind);
        X_test = Fless_mat_pruned(ind,:);
        y_test = oid(ind);
        mdl = svmtrain(y_train, X_train, ['-s 1 -t 2 -q -b 1']);
        [a,p,b] = svmpredict( y_test,  X_test, mdl,['-q -b 1']);
%         length(unique(a))
%         length(unique(y_test))
        corrmod(folderr,ss) = p(1);
        var = zeros(1,length(a));
        for zz=1:length(a)
            sdiff = setdiff(oid_,a(zz));
            var(zz) = b(zz,a(zz))-mean(b(zz,sdiff));
        end
        corrmod2(folderr,ss) = mean(var);
    end
end
end

figure('Position',[0 0 320 240])
hold on
bar(mean(corrmod2,1))
% yline(100/160)
xticks(1:3)
xticklabels({'S1','S2','S3'})
ylabel('Performance (%)')
savefig('svm')
print('svm','-dpng')
clear unity fless_mat fless_mat_unn Fless_mat utl_mask
save('svmpred')

%% Basic decoding - DTW
corrmat_dtw = true;
if corrmat_dtw
corrmoda = [];
corrmodb = [];
corrmodp = [];
for ss = 1:length(dirs)
    load(fullfile(dirs{ss},'sfp_feats.mat'))
    Fless_mat = vertcat(fless_mat{:});
  % Fless_mat = vertcat(feat_mat{:});
    anatdir = fullfile('D:\Work\ARC\',sprintf('ARC%02d',ss),'single');
    
    if ss==3; s2 = 4; else; s2 = ss; end
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
    
    Fless_mat_pruned = Fless_mat(:,1:75);
    % Fless_mat_pruned = Fless_mat(:,[3 4 9:14]);
    % [p_val, actual_statistic, permuted_stats,fig] = SFP_computeDTW_perm(T1, T2, optout, n_perms, figs)

    DTWout = SFP_pattern_DTW(Fless_mat_pruned , group_vec,350);
    DTWmean_ = cellfun(@mean,DTWout,'UniformOutput',false);
    DTWmean_ = vertcat(DTWmean_{:});
    DTWmean(:,1) = (DTWmean_(:,1)+DTWmean_(:,2))/2;
    DTWmean(:,2) = DTWmean_(:,3);
    corrmoda(:,:,ss) = DTWmean;
    % corrmodp(ss) = max(signrank(DTWmean(:,1),DTWmean(:,3)),signrank(DTWmean(:,2),DTWmean(:,3)));
end
end

figure('Position',[0.5 0.5 320 240])
bars = squeeze(mean(mean(corrmoda,3),1));
std_ = squeeze(std(mean(corrmoda,1),[],3))./sqrt(3);
bar(bars)
hold on
errorbar(1:2,bars,std_,'.')
    vec = squeeze(mean(corrmoda,1));
for ss = 1:3

    plot([1:2],vec(:,ss));
end