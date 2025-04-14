%% RSA time-series

dwnsample = 100;
switcher = 'basic'; %'basic', 'int' or 'pls';

nodor = 160;
wind = 75; % Number of samples
t_idx = 1:dwnsample:7500;
dirs = {'C:\Work\SFP\sfp_behav_s01_correct';
    'C:\Work\SFP\sfp_behav_s02_correct';
    'C:\Work\SFP\sfp_behav_s04_correct'};
%   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));
ndesc = 18;

hold on
corrmod = zeros(3,wind,ndesc );
corrmod_t = zeros(3,wind,ndesc );

sniff_trace = zeros(3,wind);
for ss = 1:length(dirs)
    load(fullfile(dirs{ss},'sfp_feats_main.mat'))
    load(fullfile(dirs{ss},'sfp_feats_main.mat'))
    load(fullfile(dirs{ss},'task_struct_trialwise.mat'))

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
    % Fless_mat_pruned = zscore(Fless_mat_pruned,[],1);

    % Behavioral features
    behav_ratings = behav.behav(ss).ratings;
    behav_ratings_ =  behav_ratings(group_vec,:);
    behav_ratings_  = zscore(  behav_ratings_ ,[],1);
    
    win = 50;
    for pp = 1:ndesc
        b1 = -(abs(behav_ratings_(:,pp)-behav_ratings_(:,pp)')');
        % b1 = pdist(behav_ratings_);
        % b1 = squareform(b1);
        DM_mat = [b1(utl_mask) unity(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)];
        % DM_mat2 = [b1(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)];

        fprintf('Subject: %02d, desc: %02d',ss, pp)
        fprintf('\n')
        for jj=1:wind; fprintf('.'); end
        fprintf('\n')
        for jj = 1:wind
            fprintf('|')
            % Fless_corr = -(abs(Fless_mat_pruned(:,jj)-Fless_mat_pruned(:,jj)'));

            % Fless_corr = corrcoef(Fless_mat(:,t_idx(jj):t_idx(jj)+win)');
            Fless_corr = -pdist(Fless_mat(:,t_idx(jj):t_idx(jj)+win));
            Fless_corr = squareform(Fless_corr);

            [wt2,t_sc2] = ARC_multicomputeWeights_tsc(DM_mat, Fless_corr(utl_mask));
            corrmod(ss,jj,pp) = wt2(2);
            corrmod_t(ss,jj,pp) = t_sc2(2);
        end
    end
    sniff_trace(ss,:) = mean(Fless_mat_pruned ,1);
    fprintf('\n')
end

%
% figure('Position',[0.5 0.5 1280 320])

SFP_clearLargeVariables
save('multidisc2')

% Figure
figure()
hold on
for ss = 1:3
    taxis = (0:wind-1)/10;
subplot(1,3,ss)
hold on
tempmat = squeeze(corrmod(ss,:,:))';
% [~,argsortmax] = max(tempmat');
% [~,argsort] = sort(argsortmax); 
argsort = 18:-1:1;
imagesc(taxis,1:18,tempmat(argsort,:))
yticks(1:18)
yticklabels(behav.behav(ss).percepts(argsort))
axis tight
colorbar
clim([-0.03 0.08])
end
savefig('disc')