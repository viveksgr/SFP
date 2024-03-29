%% Trialwise nuisance regression
root = 'C:\Work';
    dirs2 = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
        fullfile(root ,'\SFP\sfp_behav_s02_correct');
        fullfile(root ,'\SFP\sfp_behav_s04_correct')};

s = 1;
if s==1; runvec = [100*ones(1,24) 90*ones(1,24)]; else; runvec = 90*ones(1,48); end
task_run = [];
for tt = 1:length(runvec)
    task_run = blkdiag(task_run,ones(runvec(tt)));
end

if s==1; sessvec = [400*ones(1,6) 360*ones(1,6)]; else; sessvec = 360*ones(1,12); end
sess_run = [];
for tt = 1:length(sessvec)
    sess_run = blkdiag(sess_run,ones(sessvec(tt)));
end

if s==1; setvec = [1200*ones(1,2) 3*360*ones(1,2)]; else; setvec = 3*360*ones(1,4); end
set_run = [];
for tt = 1:length(setvec)
    set_run = blkdiag(set_run,ones(setvec(tt)));
end

save(fullfile(dirs2{s},'task_struct_trialwise.mat'),'task_run','sess_run','set_run')