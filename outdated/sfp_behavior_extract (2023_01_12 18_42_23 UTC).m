% Create SPM.mat for odor detection - check activity in pyriform
s = 4;
n_volumes = 876; % Number of volumes 856 for S1, 876 for S2, S3;
sn = sprintf('NEMO_%02d', s);
set_i = 1; % Initial set to be used in training
set_f = 4;  % Final set
sess_i = 2; % Initial session (for all sets)
sess_f = 4; % Final session (for all sets)
% Change filepaths here:
modelname = 'sfp_behav_s03';
create_sfpfeats = false; % isolate sniffing features from raw data
grp_ = false; % Group level analyses
corrmat_ = true;
fless_maxsize = 70000;
dwnsample = 100;

 
 
rootf = 'C:\Data\SFP';
root = 'C:\Data\NEMO';
addpath(genpath('C:\Toolboxes\breathmetrics-master'))
statpath = fullfile(rootf, modelname);
behav = load(fullfile('C:\Data\ARC\ARC','NEMO_perceptual.mat'));

% Extract valence info:
% behav_ratings = behav.behav(s).ratings(:,2);
% behav_ratings = normalize(behav_ratings,'medianiqr');

runs = 1:4; % Number of runs/session
nruns = 4*(sess_f-sess_i+1)*(set_f-set_i+1);

%% Load/Create SPM.mat
if create_sfpfeats 
mkdir(statpath)
save(fullfile(statpath,'settings.mat'));
r = 0;
onsets = cell(1,nruns);
feat_mat = cell(1,nruns);
fless_mat = cell(1,nruns);
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
    

            t2 = bmObj.time(bmObj.inhaleOnsets);           
            t_star = find_nearest(t2,onset_set);
            t_on = bmObj.inhaleOnsets(t_star);
            t_off = bmObj.exhaleOffsets(t_star);
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
            feat_mat{r} =  feat_superset(:,t_star)';
            
%           computing the sniff information - featureless extrapolated, only useful for
%           plotting
            fless_superset = zeros(1,fless_maxsize);
            for trial = 1:length(onset_set)
                n_diff = abs((t_off(trial)-t_on(trial)+1)-size(fless_superset,2));
                if size(fless_superset,2)>=(t_off(trial)-t_on(trial)+1)
                    fless_nexttrial = cat(2,R(t_on(trial):t_off(trial))',R(t_off(trial))*ones(1,n_diff));
                    fless_superset = cat(1,fless_superset,fless_nexttrial);
                else
                    fless_nexttrial = R(t_on(trial):t_off(trial))';
                    fless_prevtrial = cat(2,fless_superset,fless_superset(:,end).*ones(size(fless_superset,1),n_diff));
                    fless_superset = cat(1,fless_prevtrial ,fless_nexttrial);
                end
            end
            fless_superset(1,:) = [];
            fless_mat{r} = fless_superset(:,1:dwnsample:end);
            
%           computing the sniff information - featureless 
            fless_superset = zeros(1,fless_maxsize);
            for trial = 1:length(onset_set)
                n_diff = abs((t_off(trial)-t_on(trial)+1)-size(fless_superset,2));
                if size(fless_superset,2)>=(t_off(trial)-t_on(trial)+1)
                    fless_nexttrial = cat(2,R(t_on(trial):t_off(trial))',nan(1,n_diff));
                    fless_superset = cat(1,fless_superset,fless_nexttrial);
                else
                    fless_nexttrial = R(t_on(trial):t_off(trial))';
                    fless_prevtrial = cat(2,fless_superset,nan(size(fless_superset,1),n_diff));
                    fless_superset = cat(1,fless_prevtrial ,fless_nexttrial);
                end
            end
            fless_superset(1,:) = [];
            fless_mat_unn{r} = fless_superset(:,1:dwnsample:end);          
        end
    end
end

[max_breaths] = cellfun(@(x) size(x,2),fless_mat_unn);
fless_maxsize = (max(max_breaths)+1)*dwnsample;
fprintf('Fless Maxsize: %06d',fless_maxsize);
save(fullfile(statpath,'sfp_feats.mat'),'feat_mat','fless_mat','fless_mat_unn','onsets')
end

%% Correlation with behavior
% Sniffing modulation
if grp_
    nodor = 160;
    wind = 75; % Number of samples
    dirs = {'C:\Data\SFP\sfp_behav_s01';
            'C:\Data\SFP\sfp_behav_s02';
            'C:\Data\SFP\sfp_behav_s03'};
    color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
    behav = load(fullfile('C:\Data\ARC\ARC','NEMO_perceptual2.mat'));
    behav_ids = [1 1 1; 2 2 2; 14 13 13; 10 10 10; 6 6 6; 0 18 18];
    behav_names = behav.behav(2).percepts(behav_ids(:,2));
    
    
    figure('Position',[0 0 1280 720])
    hold on
    plot_id = 0;
    for ss = 1:length(dirs)
         load(fullfile(dirs{ss},'sfp_feats.mat'))
        Fless_mat = vertcat(fless_mat{:});
       
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
        
        for pp = 1:size(behav_ids,1)
            plot_id = plot_id +1;
            v_id = behav_ids(pp,ss);
            if v_id>0
            behav_ratings = behav.behav(ss).ratings(:,v_id);
            behav_ratings = normalize(behav_ratings,'medianiqr');
            behav_ratings_= behav_ratings(group_vec);
            divs = quantile(behav_ratings_,[0.1 0.5 0.9]);
            
            subplot(3,size(behav_ids,1),plot_id)
            hold on
            xt = 0:0.1:(wind/10-0.1);
            v_neg_hard = Fless_mat(behav_ratings_<= divs(1),1:wind );
            v_neg_soft = Fless_mat(behav_ratings_< divs(2),1:wind );
            v_pos_soft = Fless_mat(behav_ratings_>= divs(2),1:wind );
            v_pos_hard = Fless_mat(behav_ratings_>= divs(3),1:wind );
            shaded_plot(xt,mean(v_neg_hard,1),1.96*std(v_neg_hard,1)/sqrt(size(v_neg_hard,1)),color(1,:));
%             shaded_plot(xt,mean(v_neg_soft,1),1.96*std(v_neg_soft,1)/sqrt(size(v_neg_soft,1)),color(2,:));
%             shaded_plot(xt,mean(v_pos_soft,1),1.96*std(v_pos_soft,1)/sqrt(size(v_pos_soft,1)),color(3,:));
            shaded_plot(xt,mean(v_pos_hard,1),1.96*std(v_pos_hard,1)/sqrt(size(v_pos_hard,1)),color(4,:));
            legend({'low','high'})
            title(sprintf('S: %02d, %s',ss,behav_names{pp}))
            end
            if and(s==3,pp==3); xlabel('Time(s)'); end          
            if and(s==2,pp==1); ylabel('Amplitude'); end               
        end
    end   
    savefig(fullfile(rootf,'snifftraces'))
    print(fullfile(rootf,'snifftraces'),'-dpng')
end

%% Correlations Fless
if corrmat_
    nodor = 160;
    wind = 75; % Number of samples
    dirs = {'C:\Data\SFP\sfp_behav_s01';
            'C:\Data\SFP\sfp_behav_s02';
            'C:\Data\SFP\sfp_behav_s03'};
%   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
    behav = load(fullfile('C:\Data\ARC\ARC','NEMO_perceptual2.mat'));

    hold on
    corrmod = zeros(3,5);
    for ss = 1:length(dirs)
        load(fullfile(dirs{ss},'sfp_feats.mat'))
        Fless_mat = vertcat(fless_mat{:});       
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
        
        Fless_mat_pruned = Fless_mat(:,1:wind);
        Fless_corr = corrcoef(Fless_mat_pruned');
        
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

%% Correlations Fmat
 corrmat_ = true;
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
