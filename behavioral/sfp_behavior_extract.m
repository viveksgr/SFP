% Extracting sniff data from subjects and parametrizing the sniff traces
% along sniff features in each subject
% Pre-requisities:
% 1. Raw files stored in root
% 2. Behavioral data extracted in variable <behav>

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
                % if s==1
                %     R = [0; diff(R)];
                % end
                bmObj = breathmetrics(R, samplerate(1), 'humanAirflow');
                bmObj.estimateAllFeatures();

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

%% Raw Sniff plots
raw_sniff_plot = false;
if raw_sniff_plot
% Use ss = 2 and extract group_vec, Fless_mat_pruned for s2
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
end

%% Illustrations of sniffing modulation
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
        load(fullfile(dirs{ss},'sfp_feats_main.mat'))
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

