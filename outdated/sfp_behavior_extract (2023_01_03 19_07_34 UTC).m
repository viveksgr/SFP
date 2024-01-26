% Create SPM.mat for odor detection - check activity in pyriform
s = 1;
n_volumes = 856; % Number of volumes 856 for S1, 876 for S2, S3;
sn = sprintf('NEMO_%02d', s);
set_i = 1; % Initial set to be used in training
set_f = 4;  % Final set
sess_i = 1; % Initial session (for all sets)
sess_f = 3; % Final session (for all sets)

% Change filepaths here:
modelname = 'sfp_behav_s01';
rootf = 'C:\Data\SFP\Data';
root = 'C:\Data\NEMO';
addpath(genpath('C:\Toolboxes\breathmetrics-master'))
statpath = fullfile(rootf, sn, 'imaging', '1stlevelmodels', modelname);
behav = load(fullfile('C:\Data\ARC\ARC','NEMO_perceptual.mat'));

% Extract valence info:
behav_ratings = behav.behav(s).ratings(:,2);
behav_ratings = normalize(behav_ratings,'medianiqr');

runs = 1:4; % Number of runs/session
nruns = 4*(sess_f-sess_i+1)*(set_f-set_i+1);

%% Load/Create SPM.mat

mkdir(statpath)
r = 0;
detect = cell(1,nruns);
onset_set = cell(1,nruns);
odor_id = cell(1,nruns);
filename = cell(n_volumes,nruns);
noise_files = cell(nruns,1);

for set_ = set_i:set_f
    for sess = sess_i:sess_f
        for rr=runs
            r= r+1;
            if s==1
                behavpath = fullfile(root, sprintf('NEMO_%02d',s), 'behavior', sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
                load(fullfile(behavpath, sprintf('NEMO_%02d_set_%02d_sess_%02d_run_%02d.mat', s, set_, sess, rr)));
                breathpath = fullfile(root, sprintf('NEMO_%02d',s), 'breathing', sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
                load(fullfile(breathpath, sprintf('NEMO_%02d_set_%02d_sess_%02d_run_%02d.mat',s, set_, sess, rr)))
                datapath = fullfile(root, sn, 'imaging', 'nii', sprintf('set_%02d', set_), sprintf('sess_%02d', sess), sprintf('run_%02d',rr));
                ns = dir(fullfile(datapath,sprintf('nusiance_regresssors_NEMO_%02d_set_%02d_sess_%02d_run_%02d.txt',s,set_,sess,rr)));
                d = res.select_odor{1};
            else
                behavpath = fullfile(root, sprintf('NEMO_%02d',s), 'behavior', 'imaging_task',sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
                load(fullfile(behavpath, sprintf('NEMO_%02d_set_%02d_sess_%02d_run_%02d.mat', s, set_, sess, rr)));
                breathpath = fullfile(root, sprintf('NEMO_%02d',s), 'breathing', 'imaging_task', sprintf('set_%02d', set_), sprintf('sess_%02d', sess));
                load(fullfile(breathpath, sprintf('NEMO%02d_set_%02d_sess_%02d_run_%02d.mat',s, set_, sess, rr)))
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
            onset_set{r} = d(:,4)/1000 + adjust; % Condition onset
            
            % breathing stuff
            ben = data(datastart(1):dataend(1));     % when the first blip appeared. The first variable in labchart
            trigger = data(datastart(5):dataend(5));  % The 5th variable in labchart. The one with multiple blips for triggers
            
            % determine first trigger
            tr = find(round(trigger./max(trigger))==1);
            tr = tr(1);
            
            % grab breathing data
            R = data(datastart(4):dataend(4));   % respiration variable from labchart,4th variable in labchart
%             
%           % take task data
            end_time = tr+ floor(n_volumes*1.4*samplerate(1));
            if length(R)>=end_time
                R = R(tr:end_time);
            else
                R_y = length(R)+1:end_time;
                R_val = R(end)*ones(1,length(R_y));
                R = [R(tr:end) R_val];
            end
%             
%           % smooth
            R = smoothdata(R,'movmean',250)';
            if s==1;
                R = [0; diff(R)];
            end
            bmObj = breathmetrics(R, samplerate(1), 'humanAirflow');
            bmObj.estimateAllFeatures();
%             fig = bmObj.plotCompositions('raw');
            fig = bmObj.plotFeatures({'extrema','maxflow'});
%             fig = bmObj.plotCompositions('normalized');
%             fig = bmObj.plotCompositions('line');
    
            t1 = bmObj.inhaleOnsets;
            t2 = bmObj.time(bmObj.inhaleOnsets);
            
            t_star = find_nearest(t2,onset_set{r});
        end
    end
end

%% Vertical concatenation of all runs
file_list = cat(1,filename(:)); % All *.nii images to be used in training

r_col = (n_volumes*1.4)*(0:1:length(onset_set)-1); % Adjust onset times for vertical concatenation
onset_times = cell(size(onset_set));
for ii = 1:length(r_col)
    onset_times{ii} = onset_set{ii}+r_col(ii);
end
onset_set = vertcat(onset_times{:});


