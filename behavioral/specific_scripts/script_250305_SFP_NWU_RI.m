%% Load in PowerLab data for respiratory and DAQ signals

% Subject RI
dat = ('C:\Work\SFP\SFP_common\BehavioralFollowUp\PowerLabData\250305_SFP_NWU_RI.mat'); % Load PowerLab raw data as .mat file
behav_dir = 'C:\Work\SFP\SFP_common\BehavioralFollowUp\behavioralData\odor_sniffMod_results_250305_SFP_NWU_RI';
run_id = [2 4 5 6 7];
onsets = load('C:\Work\SFP\SFP_common\BehavioralFollowUp\allsniffonsets_RI.mat');
onsets = onsets.allsniffonsets_RI;
k_id = 1;

% % Subject JN
dat = ('C:\Work\SFP\SFP_common\BehavioralFollowUp\PowerLabData\250305_SFP_NWU_JN.mat'); % Load PowerLab raw data as .mat file
behav_dir = 'C:\Work\SFP\SFP_common\BehavioralFollowUp\behavioralData\odor_sniffMod_results_250305_SFP_NWU_JN';
run_id = [2 3 4 5 6];
onsets = load('C:\Work\SFP\SFP_common\BehavioralFollowUp\allsniffonsets_JN.mat');
onsets = onsets.allsniffonsets_JN;
k_id = 2;

% Subject PP
dat = ('C:\Work\SFP\SFP_common\BehavioralFollowUp\PowerLabData\250312_SFP_NWU_PP.mat'); % Load PowerLab raw data as .mat file
behav_dir = 'C:\Work\SFP\SFP_common\BehavioralFollowUp\behavioralData\odor_sniffMod_results_250312_SFP_NWU_PP';
run_id = [2 3 4 5 6];
onsets = load('C:\Work\SFP\SFP_common\BehavioralFollowUp\allsniffonsets_PP.mat');
onsets = onsets.allsniffonsets_PP;
k_id = 3;

% Subject HRM
dat = ('C:\Work\SFP\SFP_common\BehavioralFollowUp\PowerLabData\250324_SFP_NWU_HRM.mat'); % Load PowerLab raw data as .mat file
behav_dir = 'C:\Work\SFP\SFP_common\BehavioralFollowUp\behavioralData\odor_sniffMod_results_250324_SFP_NWU_HRM';
run_id = [2 3 4 5 6];
onsets = load('C:\Work\SFP\SFP_common\BehavioralFollowUp\allsniffonsets_HRM.mat');
onsets = onsets.allsniffonsets_HRM;
k_id = 4;

srate = 1000; % sampling rate is 1000 Hz
[RI_data] = ReadLabChartMat(dat); % Read in LabChart mat file
for rec = 1:length(RI_data.data) % Sort out respiratory and daq channels for each recording
    RI_resp{rec} = RI_data.data{rec,1};
    RI_daq{rec} = RI_data.data{rec,2};
end

R = cat(2,RI_resp{run_id});
RI_daq_allRuns = cat(2,RI_daq{run_id});
% 
% figure()
% yyaxis  left
% plot(RI_resp_allRuns) % Plot combined respiratory signal
% yyaxis  right
% plot(RI_daq_allRuns) % Plot combined daq signal

%% Breathing data extraction
bmObj = breathmetrics(R, 1000, 'humanAirflow');
bmObj.estimateAllFeatures();

t_star = find_nearest(bmObj.inhaleOnsets,onsets')'; % Done in samples
t_on = bmObj.inhaleOnsets(t_star);

% computing the sniff information - feature wise
prop_list = properties(bmObj);
proper_list = prop_list(7:24);
feat_superset = [];
for feat = 1:length(proper_list)
    eval(sprintf('feat_superset = cat(1,feat_superset,bmObj.%s);', ...
        proper_list{feat}))
end
feat_superset =  feat_superset(:,t_star)';

% peak times
peak_times = bmObj.inhalePeaks(t_star);
assert(sum(peak_times<t_on)==0,'Peak times must always be less than inhale onsets' )
fless_size = 7500;
fless_superset = zeros(length(onsets),  fless_size);
for trial = 1:length(onsets)
    % Use this for odor-onset aligned
    % fless_superset(trial,:) = R( R_times_idx(trial):R_times_idx(trial)+fless_size-1)';

    % Use this for sniffonset aligned
    fless_superset(trial,:) = R( t_on(trial):t_on(trial)+fless_size-1)';
end
fless_mat{k_id} = fless_superset;
features = SFP_extractSniffFeatures(fless_superset, feat_superset,1000);
feat_mat{k_id} =  cat(2,feat_superset(),features);


%% Load in behavioral data
behav_dir_list = dir(fullfile(behav_dir,'*.mat'));

if k_id~=3
assert(length(behav_dir_list)==5) % Only 5 runs
odor_order = [];
for zz=1:5 
    load(fullfile(behav_dir,behav_dir_list(zz).name)); 
    odor_order = [odor_order task_vars.odor_order(1:40)]; 
end
else
    temp_mat = load(fullfile(behav_dir,'taskVars_250312_SFP_NWU_PP.mat'));
    oid_runs = temp_mat.taskVars_250305_SFP_NWU_PP;
    odor_order = [];
    for zz=1:5
        eval(sprintf('temp_runs=oid_runs.run%01d.odor_order;',zz))
        odor_order = [odor_order temp_runs];
    end
end
odor_id{k_id} = odor_order';
save('sfp_data_temp.mat','feat_mat','fless_mat','odor_id')

% 
% load('/Volumes/Lab_Common/SFP_common/BehavioralFollowUp/behavioralData/odor_sniffMod_results_250305_SFP_NWU_RI/250305_SFP_NWU_RI_results_run1.mat','results','task_vars'); % load run 1 results
% results_run1 = results;
% task_vars_run1 = task_vars;
% load('/Volumes/Lab_Common/SFP_common/BehavioralFollowUp/behavioralData/odor_sniffMod_results_250305_SFP_NWU_RI/250305_SFP_NWU_RI_results_run2.mat','results','task_vars'); % load run 2 results
% results_run2 = results;
% task_vars_run2 = task_vars;
% load('/Volumes/Lab_Common/SFP_common/BehavioralFollowUp/behavioralData/odor_sniffMod_results_250305_SFP_NWU_RI/250305_SFP_NWU_RI_results_run3.mat','results','task_vars'); % load run 3 results
% results_run3 = results;
% task_vars_run3 = task_vars;
% load('/Volumes/Lab_Common/SFP_common/BehavioralFollowUp/behavioralData/odor_sniffMod_results_250305_SFP_NWU_RI/250305_SFP_NWU_RI_results_run4.mat','results','task_vars'); % load run 4 results
% results_run4 = results;
% task_vars_run4 = task_vars;
% load('/Volumes/Lab_Common/SFP_common/BehavioralFollowUp/behavioralData/odor_sniffMod_results_250305_SFP_NWU_RI/250305_SFP_NWU_RI_results_run5.mat','results','task_vars'); % load run 5 results
% results_run5 = results;
% task_vars_run5 = task_vars;
% 
% %% Combine ratings, odor odor, 
% ratings_allRuns = cat(2,results_run1.rating,results_run2.rating,results_run3.rating,...
%     results_run4.rating,results_run5.rating);
% rating_RT_allRuns = cat(2,results_run1.ratingRT,results_run2.ratingRT,results_run3.ratingRT,...
%     results_run4.ratingRT,results_run5.ratingRT)/srate;
% odor_order_allRuns = cat(2,task_vars_run1.odor_order,task_vars_run2.odor_order,...
%     task_vars_run3.odor_order(1:40),task_vars_run4.odor_order,task_vars_run5.odor_order);
% desc_order_allRuns = cat(2,task_vars_run1.desc_order,task_vars_run2.desc_order,...
%     task_vars_run3.desc_order(1:40),task_vars_run4.desc_order,task_vars_run5.desc_order);
% descriptors_runs1and2 = task_vars_run1.descriptors;
% descriptors_run3thru5 = task_vars_run3.descriptors;
% detect_rating_allRuns = cat(2,results_run1.detect_code,results_run2.detect_code,...
%     results_run3.detect_code,results_run4.detect_code,results_run5.detect_code);
% 
% CID = [240;8785;7335;7410;15510;7720;1068;9862;10364;5365027];
% 
% 
% %% Combine behavioral data matrix, similar to the one commented out below that Vivek used
% behav_data_mat(:,1) = odor_order_allRuns;
% behav_data_mat(:,2) = desc_order_allRuns;
% behav_data_mat(:,3) = detect_rating_allRuns;
% behav_data_mat(:,4) = ratings_allRuns;
% behav_data_mat(:,5) = rating_RT_allRuns;
% 
% 
% 
%     % Data array
%     % [1: odor order (1-10)
%     %  2: CID of odor
%     %  3: onset of odor trigger (sec)
%     %  4: onset of sniff cue (sec)
%     %  5: descriptor ID
%     %  6: detect rating
%     %  7: button pressed (mouse left/right)
%     %  8: detect RT (s)
%     %  9: time at which detection was submitted
%     % 10: descriptor rating
%     % 11: descriptor RT (s)
%     % 12: time at which descriptor rating was submitted
%     % 13: (optional) initial starting point of scale
%     %  -- You can expand these columns as needed.

