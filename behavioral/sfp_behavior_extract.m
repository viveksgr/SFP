% Create SPM.mat for odor detection - check activity in pyriform
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
                %
                %           % smooth
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



                % % % % Old method
                % % % t2 = bmObj.time(bmObj.inhaleOnsets);
                % % % t_star = find_nearest(t2,onset_set);
                % % % t_on2 = bmObj.inhaleOnsets(t_star);
                % % % if adjuster
                % % %     load(fullfile(statpath,'corr_adjust',sprintf('Adjust_onsets_Subject_%01d-Set_%01d-Session_%01d-Run_%01d.mat',s,set_,sess,rr)),'adjust_onsets')
                % % %     t_on = adjust_onsets';
                % % % else
                % % %     t_on = t_on2;
                % % % end
                
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

                % 
                % % % computing the sniff information - featureless extrapolated, only useful for
                % % % plotting
                % fless_superset = zeros(1,fless_maxsize);
                % for trial = 1:length(onset_set)
                %     n_diff = abs((t_off(trial)-t_on(trial)+1)-fless_maxsize);
                %     if fless_maxsize>=(t_off(trial)-t_on(trial)+1)
                %         fless_nexttrial = cat(2,R(t_on(trial):t_off(trial))',R  (t_off(trial))*ones(1,n_diff));                   
                %     else
                %         fless_nexttrial = R(t_on(trial):t_on(trial)+fless_maxsize-1)';
                %     end
                %      fless_superset = cat(1,fless_superset,fless_nexttrial);
                % end
                % fless_superset(1,:) = [];
                % fless_mat{r} = fless_superset(:,1:dwnsample:end);


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
grp_ = false;
rootf = 'C:\Work\SFP\fig_plots\behavioral_panels';
if grp_
    nodor = 160;
    wind = 75; % Number of samples
    dirs = {'C:\Work\SFP\sfp_behav_s01';
        'C:\Work\SFP\sfp_behav_s02';
        'C:\Work\SFP\sfp_behav_s03'};
    color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
    behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));
    behav_ids = [1 1 1; 2 2 2; 14 13 13; 10 10 10; 6 6 6; 3 3 3; 4 4 4; 5 5 5];
    behav_names = behav.behav(2).percepts(behav_ids(:,2));


    figure('Position',[0 0 1280 720])
    hold on
    plot_id = 0;
    for ss = 1:length(dirs)
        load(fullfile(dirs{ss},'sfp_feats.mat'))
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
                xt = 0:0.1:(wind/10-0.1);
                % v_neg_hard = Fless_mat(behav_ratings_<= divs(1),1:wind );
                % v_pos_hard = Fless_mat(behav_ratings_>= divs(3),1:wind );

                v_neg_hard = Fless_mat(behav_ratings_<= divs2(2),1:wind );
                v_pos_hard = Fless_mat(behav_ratings_>= divs2(end-1),1:wind );
                v_null_hard = Fless_mat(and(behav_ratings_>= divs2(9),behav_ratings_<= divs2(11)),1:wind );


                feat_neg_hard = feat_mat(behav_ratings_<= divs(1),[3 4 9:14]);
                feat_pos_hard = feat_mat(behav_ratings_>= divs(3),[3 4 9:14]);
                % 
                % shaded_plot(xt,mean(v_neg_hard,1),1.96*std(v_neg_hard,1)/sqrt(size(v_neg_hard,1)),color(1,:));
                % %             shaded_plot(xt,mean(v_neg_soft,1),1.96*std(v_neg_soft,1)/sqrt(size(v_neg_soft,1)),color(2,:));
                % %             shaded_plot(xt,mean(v_pos_soft,1),1.96*std(v_pos_soft,1)/sqrt(size(v_pos_soft,1)),color(3,:));
                % shaded_plot(xt,mean(v_null_hard,1),1.96*std(v_null_hard,1)/sqrt(size(v_null_hard,1)),color(2,:));
                % shaded_plot(xt,mean(v_pos_hard,1),1.96*std(v_pos_hard,1)/sqrt(size(v_pos_hard,1)),color(4,:));

                imagesc(Fless_mat(argsort,1:wind))
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
                    legend({'low','mid','high'})
                end
                title(sprintf('S: %02d, %s',ss,behav_names{pp}))
            end
            if and(ss==3,pp==3); xlabel('Time(s)'); end
            if and(ss==2,pp==1); ylabel('Amplitude'); end
        end
    end
    savefig(fullfile(rootf,'snifftraces_mat'))
    print(fullfile(rootf,'snifftraces_mat'),'-dpng')
end

%% Correlations Fless
if corrmat_
    nodor = 160;
    wind = 75; % Number of samples
    dirs = {'C:\Work\SFP\sfp_behav_s01';
        'C:\Work\SFP\sfp_behav_s02';
        'C:\Work\SFP\sfp_behav_s03'};
    %   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
    behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));

    hold on
    corrmod = zeros(3,5);
    for ss = 1:length(dirs)
        load(fullfile(dirs{ss},'sfp_feats_corr.mat'))
        fless_mat{1,19}(:,701:end) = [];
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

        Fless_mat_pruned = Fless_mat(:,1:wind);
        Fless_corr = corrcoef(Fless_mat_pruned');
        Fless_corr(isnan(Fless_corr))=0;

        % Behavioral features
        behav_ratings = behav.behav(ss).ratings;
        behav_ratings_ =  behav_ratings(group_vec,:);

        b1 = corrcoef(behav_ratings_')';
        temp = corrcoef(b1,Fless_corr);
        corrmod(ss,1) = temp(2);

        bint = corrcoef(behav_ratings_(:,2:end)')';
        temp = corrcoef(bint,Fless_corr);
        corrmod(ss,2) = temp(2);

        bpls = corrcoef(behav_ratings_(:,[1 3:end])')';
        temp = corrcoef(bpls,Fless_corr);
        corrmod(ss,3) = temp(2);

        rpls = regressmeout(behav_ratings(:,[1 3:end])',repmat(behav_ratings(:,2)',17,1));
        rpls = rpls(:,group_vec);
        bpls = corrcoef(rpls);
        temp = corrcoef(bpls,Fless_corr);
        corrmod(ss,4) = temp(2);
        
        rint = regressmeout(behav_ratings(:,[2:end])',repmat(behav_ratings(:,1)',17,1));
        rint = rint(:,group_vec);
        brint = corrcoef(rint);
        temp = corrcoef(brint,Fless_corr);
        corrmod(ss,5) = temp(2);

     
    end


    figure('Position',[0.5 0.5 640 480])
    bar(mean(corrmod))
    hold on
    errorbar([1:5],mean(corrmod),std(corrmod)/sqrt(3),'.')
    cs = {'r','g','b'};
    for ss = 1:3
        plot([1:5],corrmod(ss,:),cs{ss})
    end
    yline(r2t(0.05,nchoosek(length(utl_mask),2)))
    ylabel('representational similarity')
    xticks([1:5])
    xticklabels({'All','-Int','-Pls','-Pls (reg)','- Int (reg)',})
end

%% Correlations Fmat

%  corrmat_ = true;
if corrmat_
    nodor = 160;
    dirs = {'C:\Data\SFP\sfp_behav_s01';
        'C:\Data\SFP\sfp_behav_s02';
        'C:\Data\SFP\sfp_behav_s03'};
    %   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
    behav = load(fullfile('C:\Data\ARC\ARC','NEMO_perceptual2.mat'));

    hold on
    corrmod = zeros(3,5);
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
        %                 sum(isnan( Fless_mat_pruned(:)))
        Fless_corr = pdist(Fless_mat_pruned);
        Fless_corr = 1-squareform(Fless_corr);

        %         Fless_corr = corrcoef(Fless_mat_pruned','rows','pairwise');

        % Behavioral features
        behav_ratings = behav.behav(ss).ratings;
        behav_ratings_ =  behav_ratings(group_vec,:);

        b1 = corrcoef(behav_ratings_')';
        temp = corrcoef(b1,Fless_corr);
        corrmod(ss,1) = temp(2);

        bint = corrcoef(behav_ratings_(:,2:end)')';
        temp = corrcoef(bint,Fless_corr);
        corrmod(ss,2) = temp(2);

        bpls = corrcoef(behav_ratings_(:,[1 3:end])')';
        temp = corrcoef(bpls,Fless_corr);
        corrmod(ss,3) = temp(2);

        rint = regressmeout(behav_ratings(:,[2:end])',repmat(behav_ratings(:,1)',17,1));
        rint = rint(:,group_vec);
        brint = corrcoef(rint);
        temp = corrcoef(brint,Fless_corr);
        corrmod(ss,4) = temp(2);

        rpls = regressmeout(behav_ratings(:,[1 3:end])',repmat(behav_ratings(:,2)',17,1));
        rpls = rpls(:,group_vec);
        bpls = corrcoef(rpls);
        temp = corrcoef(bpls,Fless_corr);
        corrmod(ss,5) = temp(2);
    end
end

figure()
bar(mean(corrmod))
hold on
errorbar([1:5],mean(corrmod),std(corrmod)/sqrt(3),'.')
cs = {'r','g','b'};
for ss = 1:3
    plot([1:5],corrmod(ss,:),cs{ss})
end
yline(r2t(0.05,nchoosek(length(utl_mask),2)))
ylabel('representational similarity')
xticks([1:5])
xticklabels({'All','-Int','-Pls','Int reg','Pls reg'})
savefig('Feats')
print('Feats','-dpng')

%% RSA time-series
corrmat_ = true;
dwnsample = 100;
if corrmat_
    nodor = 160;
    wind = 75; % Number of samples
    dirs = {'C:\Work\SFP\sfp_behav_s01_correct';
        'C:\Work\SFP\sfp_behav_s02_correct';
        'C:\Work\SFP\sfp_behav_s03_correct'};
    %   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
    behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));

    hold on
    corrmod = zeros(3,wind);
    corrmod_int = zeros(3,wind);
    corrmod_pls = zeros(3,wind);
    sniff_trace = zeros(3,wind);
    for ss = 1:length(dirs)
        ss
        load(fullfile(dirs{ss},'sfp_feats_corrected.mat'))
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