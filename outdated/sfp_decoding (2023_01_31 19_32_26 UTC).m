%% Find raw decodability

nodor = 160;
wind = 75; % Number of samples
dirs = {'C:\Data\SFP\sfp_behav_s01';
    'C:\Data\SFP\sfp_behav_s02';
    'C:\Data\SFP\sfp_behav_s03'};
%   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
behav = load(fullfile('C:\Data\ARC\ARC','NEMO_perceptual2.mat'));

hold on
if corrmat_
corrmod = zeros(3,2);
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
    
    
    M_on = logical(unity);
    M_on(logical(eye(size(M_on)))) = false;
    
    M_off = ~logical(unity);
    % Behavioral features
    
    pval = ranksum(Fless_corr(M_on),Fless_corr(M_off));    
    corrmod(ss,1) = mean(Fless_corr(M_on))-mean(Fless_corr(M_off));   
    corrmod(ss,2) =  pval;
end
end
%% Make SVM





%% Compute if errors are for perceptually similar odors