s = 4;
n_volumes = 876; % Number of volumes 856 for S1, 876 for S2, S3;
sn = sprintf('NEMO_%02d', s);
set_i = 1; % Initial set to be used in training
set_f = 4;  % Final set
sess_i = 2; % Initial session (for all sets)
sess_f = 4; % Final session (for all sets)

% Change filepaths here:
modelname = 'sfp_behav_s04_correct';
create_sfpfeats = true; % isolate sniffing features from raw data
grp_ = false; % Group level analyses
corrmat_ = false;
fless_maxsize = 7500;
dwnsample = 100;
adjuster = true;

rootf = 'C:\Work\SFP';
root = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Neurology\Kahnt_Lab\Vivek\KData\NEMO';
% addpath(genpath('C:\Toolboxes\breathmetrics-master'))
statpath = fullfile(rootf, modelname);
behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual.mat'));

% Extract valence info:
% behav_ratings = behav.behav(s).ratings(:,2);
% behav_ratings = normalize(behav_ratings,'medianiqr');

runs = 1:4; % Number of runs/session
nruns = 4*(sess_f-sess_i+1)*(set_f-set_i+1);
sniff_c = 0;
sniff_c2 = 0;
sniff_c3 = 0;

sniff_cmat = logical([]);
% Load/Create SPM.mat
if create_sfpfeats
    mkdir(statpath)
    save(fullfile(statpath,'settings.mat'));
    r = 0;
    onsets = cell(1,nruns);
    feat_mat = cell(1,nruns);
    fless_mat = cell(1,nruns);
    fless_mat_peak = cell(1,nruns);
    fless_mat_unn = cell(1,nruns);
    breath_id = [4 6 0 4]; % Which channel in the labchart is breathing trace
    for set_ = set_i:set_f
        for sess = sess_i:sess_f
            for rr=runs
                r= r+1;
                if s==1
                    behavpath = fullfile(root, sprintf('NEMO_%02d',s), 'behavior', sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
                    load(fullfile(behavpath, sprintf('NEMO_%02d_set_%02d_sess_%02d_run_%02d.mat', s, set_, sess, rr)));
                    breathpath = fullfile(root, sprintf('NEMO_%02d',s), 'breathing', sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
                    load(fullfile(breathpath, sprintf('NEMO_%02d_set_%02d_sess_%02d_run_%02d.mat',s, set_, sess, rr)))
                    load(fullfile(breathpath, sprintf('time_adjust_NEMO_%02d_set_%02d_sess_%02d_run_%02d.mat',s, set_, sess, rr)))
                    datapath = fullfile(root, sn, 'imaging', 'nii', sprintf('set_%02d', set_), sprintf('sess_%02d', sess), sprintf('run_%02d',rr));
                    ns = dir(fullfile(datapath,sprintf('nusiance_regresssors_NEMO_%02d_set_%02d_sess_%02d_run_%02d.txt',s,set_,sess,rr)));
                    d = res.select_odor{1};
                else
                    behavpath = fullfile(root, sprintf('NEMO_%02d',s), 'behavior', 'imaging_task',sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
                    load(fullfile(behavpath, sprintf('NEMO_%02d_set_%02d_sess_%02d_run_%02d.mat', s, set_, sess, rr)));
                    breathpath = fullfile(root, sprintf('NEMO_%02d',s), 'breathing', 'imaging_task', sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
                    load(fullfile(breathpath, sprintf('NEMO%02d_set%02d_sess%02d_run%02d.mat',s, set_, sess, rr)))
                    load(fullfile(breathpath, sprintf('time_adjust_NEMO%02d_set_%02d_sess_%02d_run_%02d.mat',s, set_, sess, rr)))
                    datapath = fullfile(root, sn, 'imaging', 'nii', sprintf('set_%02d', set_), sprintf('sess_%02d', sess), sprintf('run_%02d',rr));
                    ns = dir(fullfile(datapath,sprintf('nuisance_regresssors_NEMO%02d_set_%02d_sess_%02d_run_%02d.txt',s,set_,sess,rr)));
                    d = res.select_odor;
                end
                % Col.  Description
                % 1     odor order (1-10)
                % 2     CID of odor
                % 3     onset of odor trigger from t0
                % 4     onset of sniff cue from t0
                % 5     percept to be rated
                % 6     Detect rating (0 = no smell)
                % 7     Mouse button pressed
                % 8     Detect RT
                % 9     Time at which 8 submitted from t0 
                % 10    Percept rating
                % 11    Percept RT
                % 12    Time at which 11 submitted from t0
                onset_set = d(:,4)/1000 + adjust; % Condition onset
                onsets{r} = onset_set;

                % breathing stuff
                ben = data(datastart(1):dataend(1));     % when the first blip appeared. The first variable in labchart
                trigger = data(datastart(5):dataend(5));  % The 5th variable in labchart. The one with multiple blips for triggers
                % determine first trigger
                tr = find(round(trigger./max(trigger))==1);
                tr = tr(1);
                % grab breathing data
                R = data(datastart(breath_id(s)):dataend(breath_id(s)));   % respiration variable from labchart,4th variable in labchart
                %           % take task data
                end_time = tr+ floor(n_volumes*1.4*samplerate(1));
                if length(R)>=(end_time)
                    R = R(tr:end_time);
                else
                    fprintf('Resp trace is clipped. Extrapolating...\n')
                    R_y = length(R)+1:end_time;
                    R_val = R(end)*ones(1,length(R_y));
                    R = [R(tr:end) R_val];
                end
                
                % smooth
                R = smoothdata(R,'movmean',250)';
                if s==1
                    R = [0; diff(R)];
                end
                bmObj = breathmetrics(R, samplerate(1), 'humanAirflow');
                bmObj.estimateAllFeatures();
                %             fig = bmObj.plotCompositions('raw');
                %             fig = bmObj.plotFeatures({'extrema','maxflow'});
                %             fig = bmObj.plotCompositions('normalized');
                %             fig = bmObj.plotCompositions('line');


                % Correct method 
                if adjuster
                    load(fullfile(statpath,'corr_adjust',sprintf('Adjust_onsets_Subject_%01d-Set_%01d-Session_%01d-Run_%01d.mat',s,set_,sess,rr)),'adjust_onsets')
                    t_star = find_nearest(bmObj.inhaleOnsets,adjust_onsets')'; % Done in samples                    
                else
                    t_star = find_nearest(bmObj.time(bmObj.inhaleOnsets),onset_set); % Done in absolute time                   
                end
                t_on = bmObj.inhaleOnsets(t_star);

                t_off = bmObj.exhaleOffsets(t_star);
                t_off_ind = t_off-t_on;
                % % % % t_on(t_off_ind<0) = t_on2(t_off_ind<0);
                assert(sum((t_off-t_on)<0,'all')==false) % t_off must be always greater than t_on

           
                if isnan(t_off(end))
                    t_off(isnan(t_off))= length(R);
                end
                fprintf('Number of overlapping traces %02d\n',sum(diff(t_star)==0))

                % computing the sniff information - feature wise
                prop_list = properties(bmObj);
                proper_list = prop_list(7:24);
                feat_superset = [];
                for feat = 1:length(proper_list)
                    eval(sprintf('feat_superset = cat(1,feat_superset,bmObj.%s);',proper_list{feat}))
                end
                feat_superset =  feat_superset(:,t_star)';

                % peak times
                peak_times = bmObj.inhalePeaks(t_star);
                assert(sum(peak_times<t_on)==0,'Peak times must always be less than inhale onsets' )


                % computing the sniff information - at peaks
                R_times = (0:1:length(R))/samplerate(1);
                R_times_idx = find_nearest(R_times,onset_set);

                % Count t_ons
                

                % 
                % sniff_c = sniff_c+sum((R_times(t_on2)'-onset_set)<-0.5);
                % % if (sum((R_times(t_on2)'-onset_set)<-0.5))>15
                % %     'beep'
                % % end
                % sniff_c2 = sniff_c2+sum((R_times(t_on)'-onset_set)<-0.5);

                sniff_cmat = cat(1,sniff_cmat,(R_times(t_on)'-onset_set)>-0.5);

                sniff_c3 = sniff_c3+sum((peak_times'-onset_set)<-0.5);


                fless_size = 7500; 
                fless_superset = zeros(length(onset_set),  fless_size);
                peaks = zeros(length(onset_set),1);
                troughs = zeros(length(onset_set),1);
                for trial = 1:length(onset_set)
                    % Use this for odor-onset aligned 
                     % fless_superset(trial,:) = R( R_times_idx(trial):R_times_idx(trial)+fless_size-1)'; 
                      
                     % Use this for sniffonset aligned
                     fless_superset(trial,:) = R( t_on(trial):t_on(trial)+fless_size-1)';
                end
                fless_mat{r} = fless_superset;
     
                features = SFP_extractSniffFeatures(fless_superset, feat_superset,1000);


                feat_mat{r} =  cat(2,feat_superset(),features);

                % % computing the sniff information - at peaks
                % fless_size = 5000; 
                % fless_superset = zeros(length(onset_set),  fless_size);
                % for trial = 1:length(onset_set)
                %      fless_superset(trial,:) = R(peak_times(trial)-1000:peak_times(trial)+fless_size-1001)';
                % end
                % fless_mat_peak{r} = fless_superset(:,1:dwnsample:end);

                % prop_list = properties(bmObj);
                % proper_list = prop_list(7:24);
                % feat_superset = [];
                % for feat = 1:length(proper_list)
                %     eval(sprintf('feat_superset = cat(1,feat_superset,bmObj.%s);',proper_list{feat}))
                % end
                % feat_mat{r} =  feat_superset(:,t_star)';

            end
        end
    end

    [max_breaths] = cellfun(@(x) size(x,2),fless_mat_unn);
    % fless_maxsize = (max(max_breaths)+1)*dwnsample;
    % fprintf('Fless Maxsize: %06d',fless_maxsize);
    fprintf('Sniff counts: %.2f',sniff_c/4320)
    fprintf('Sniff counts_adj: %.2f',sniff_c2/4320)
    fprintf('Sniff counts_peak: %.2f',sniff_c3/4320)
    save(fullfile(statpath,'sfp_feats_main.mat'),'feat_mat','fless_mat','fless_mat_unn','fless_mat_peak','onsets','sniff_cmat')
end

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

%% RSA - Trial-wise
numpcs = [14 11 11]; % 90% Variance
root = 'C:\Work';
corrmat_ = true;
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
    savepath = 'C:\Work\SFP\Final_plots\Behavioral\Trialwise RSA';
    hold on
    nconds = 4;
    rsa_P1 = zeros(3,nconds,2);
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
        A2_corr = pdist(mainmat,"spearman");
        A2_corr = 1-squareform(A2_corr);

        behav_corr = corrcoef(behav_ratings');
        int_corr = -abs(behav_ratings(:,1)-behav_ratings(:,1)');
        pls_corr = -abs(behav_ratings(:,2)-behav_ratings(:,2)');

        [~,t_sc] = ARC_multicomputeWeights_tsc([behav_corr(utl_mask) unity(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)], A2_corr(utl_mask) );
        rsa_P1(ss,1,1) = t_sc(2);
        rsa_P1(ss,1,2) = t_sc(3);

        [~,t_sc] = ARC_multicomputeWeights_tsc([behav_corr(utl_mask) unity(utl_mask) int_corr(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)],  A2_corr(utl_mask));
        rsa_P1(ss,2,1) = t_sc(2);
        rsa_P1(ss,2,2) = t_sc(3);

        [~,t_sc] = ARC_multicomputeWeights_tsc([behav_corr(utl_mask) unity(utl_mask) pls_corr(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)],  A2_corr(utl_mask));
        rsa_P1(ss,3,1) = t_sc(2);
        rsa_P1(ss,3,2) = t_sc(3);

        [~,t_sc] = ARC_multicomputeWeights_tsc([behav_corr(utl_mask) unity(utl_mask) int_corr(utl_mask)  pls_corr(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)],  A2_corr(utl_mask));
        rsa_P1(ss,4,1) = t_sc(2);
        rsa_P1(ss,4,2) = t_sc(3);
    end

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
    
    % Subject data points
    c_s = {'r','g','b'}; % Data dots for subjects
    for ii = 1:nconds % For bars for perceptual, chemical and combinations
        for jj = 1:3
            plot(x_m(:,ii),squeeze(rsa_P1(jj,ii,:)),c_s{jj},'handle','off')
        end
    end
    yline(tinv(0.95,sum(utl_mask(:))))
    xticks(1:4)
    xticklabels({'Sniff RSA','-Int','-Pls','-Int-Pls'})
    ylabel('Representational Similarity (t)')
    legend({'Perceptual similarity','Odor trial similarity'})
    % yline(r2t(0.05,sum(utl_mask2(:))));
    % yline(r2t(0.05,nchoosek(length( group_vec),2)));
    savefig(fullfile(savepath,'fless_map'))
    print(fullfile(savepath,'fless_map'),'-dpng')
    % clear Fmat_1_m behav_mat unity M_anat M_reg n_reg M_unstack fless_mat fless_mat_unn modelmd_ modelmd S_omat_vals utl_mask utl_mask2
    save(fullfile(savepath,'ARC_RSA_fless'),'settings_','rsa_P1')
end
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
    load(fullfile(dirs{ss},'sfp_feats_corrected.mat'))

  

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

%% RSA static - Cluster
clust_maker = true;
repcounts = 1; 

settings_.featorflessnot = true;
settings_.numClusters = 20;
dwn = 50;
if clust_maker
root = 'C:\Work';

nodor = 160;
wind = 3500; % Number of samples
dirs = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
        fullfile(root ,'\SFP\sfp_behav_s02_correct');
        fullfile(root ,'\SFP\sfp_behav_s04_correct')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
        fullfile(root,'ARC\ARC\ARC02\single');
        fullfile(root,'ARC\ARC\ARC03\single')};
% 
% dirs3 = {fullfile(root,'SFP\CCA_2comp_full\SFP1');
%         fullfile(root,'SFP\CCA_2comp_full\SFP2');
%         fullfile(root,'SFP\CCA_2comp_full\SFP3')};
savepath = 'C:\Work\SFP\Results\Clustering\behavioral_RSA';

nplt = 10;
%   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));
hold on

corrmod = zeros(3,3);
Amat = {};

sbplt  =0;
tic
figure()
idxmats = cell(3,1);
hold on
fprintf('\n')
for ss = 1:length(dirs)
    fprintf('Subject %02d\n',ss)
    statpath = dirs{ss};
    anatdir = dirs2{ss};

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

    behav_ratings = behav.behav(ss).ratings;
    utl_mask = logical(triu(ones(length(behav_ratings)),1)); % All possible odors
    behav_ratings2 = behav_ratings(group_vec,:);

    % Sniff descriptors
    load(fullfile(dirs{ss},'sfp_feats_corrected.mat'))
  
    Fless_mat = vertcat(fless_mat{:});
    Fless_mat_pruned = Fless_mat(:,1:dwn:wind);  
    Feat_mat_pruned = vertcat(feat_mat{:});
    Feat_mat_pruned =  Feat_mat_pruned(:,[3 4 9:14]) ;

    if settings_.featorflessnot; mainmat = Feat_mat_pruned; else; mainmat = Fless_mat_pruned; end
    mainmat(isnan(mainmat))=0;
    if settings_.featorflessnot
        mainmat = zscore(mainmat,[],1);
    end  

    D1 = (1-corrcoef(mainmat'))/2'; % Raw correlation value

    mainmatmax = cat(2,mainmat,behav_ratings2);
    % D2 = arrayfun(@(row) SFP_constraineddist( mainmatmax(row,:), mainmatmax), 1:size(mainmatmax, 1),'UniformOutput',false);
    % D2 = vertcat(D2{:});

    % [DO, i_sort, idx,sumd] = SFP_spectral_reorder_and_cluster(D, settings_.numClusters);
   
    % [idx,A1,sumd2] = kmeans(mainmat,settings_.numClusters,'MaxIter', 1000, 'Replicates',1000,'Distance','sqeuclidean');

    % Elbow
    % [w_val(ss,:),idx,A1] = SFP_kmeansElbowMethod(mainmat, 2:2:settings_.numClusters);
    % subplot(1,3,ss)
    % hold on
    % plot(2:2:settings_.numClusters,w_val(ss,:),'.')
    % xlabel('num clusters')
    % title(sprintf('sub%02d'),ss)
    % ylabel('Kmeans pfmance')

    % Raw Sniff Corr
    % [idx,~] = kmedoids(mainmat,settings_.numClusters,'Distance', 'Euclidean','Replicates',100);
    [~, ~, idx,~] = SFP_spectral_reorder_and_cluster(D1, settings_.numClusters);

    A1 = SFP_splitapply_mean(mainmat,idx);
    A1_corr = corrcoef(A1');
    M_cmat2 = SFP_splitapply_mean(behav_ratings2,idx);
    M_corr_sniff = corrcoef(M_cmat2');

    utl_mask2 = logical(triu(ones(length(unique(idx))),1)); % All possible odors

    A2 = splitapply(@mean,mainmat,group_vec); 
    A2_corr = corrcoef(A2'); 
    M_corr = corrcoef(behav_ratings');
  
  
    temp = corrcoef(M_corr_sniff(utl_mask2),A1_corr(utl_mask2),"Rows","Pairwise");
    [~,corrmod(ss,1)] = r2p(temp(2),sum(utl_mask2(:)));
  
    [~,corrmod(ss,2)] = r2p(fastcorr(M_corr(utl_mask),A2_corr(utl_mask)),sum(utl_mask(:)));

    for cc = 1:repcounts
        if repcounts>1
            [idx_adj,~] = kmedoids(mainmatmax,settings_.numClusters,'Distance', @SFP_constraineddist);
        else
            fprintf('Running non iterative mode\n')
            % [idx_adj,~] = kmedoids(mainmatmax,settings_.numClusters,'Distance', @SFP_constraineddist,'Replicates',100);
            D2 = arrayfun(@(row) SFP_constraineddist( mainmatmax(row,:), mainmatmax), 1:size(mainmatmax, 1),'UniformOutput',false);
            D2 = horzcat(D2{:});
            [~, ~, idx_adj,~] = SFP_spectral_reorder_and_cluster(D2, settings_.numClusters);

        end


    A1_adj = SFP_splitapply_mean(mainmat,idx_adj);
    A1_corr_adj = corrcoef(A1_adj');
    M_cmat2_adj = SFP_splitapply_mean(behav_ratings2,idx_adj);
    M_corr_sniff_adj = corrcoef(M_cmat2_adj');

    temp = corrcoef(M_corr_sniff_adj(utl_mask2),A1_corr_adj(utl_mask2),"Rows","Pairwise");
    [~,temp_tsc] = r2p(temp(2),sum(utl_mask2(:)));
    if or(temp_tsc<tinv(0.975,sum(utl_mask2(:))),repcounts==1)
        corrmod(ss,3) = temp_tsc;
        idxmats{ss} = idx_adj;
        fprintf('Model found: %02d',cc)
        break
    end

    end

end
% toc
% figure() = mean(corrmod,1);
% bar(m)
% hold on
% errorbar([1:size(corrmod,2)],m,std(corrmod),'.')
% c = {'r','g','b'};
% for ss = 1:3; plot([1:size(corrmod,2)],corrmod(ss,:),c{ss},'Marker','.','Markersize',10); end
% ylabel('Sniff RSA')
% % yline(r2t(0.05,nchoosek(4560,2)))
% xticks([1:size(corrmod,2)])
% xticklabels({'Sniff type','Odor type','Sniff type adj'})
% xtickangle(45)
% % yline(r2t(0.05,12720))
% 
% yline(1.65)
end