%--------------------------------------------------------------------------
% SCRIPT: sfp_decoding.m
%
% DESCRIPTION:
%   This script performs basic decoding analysis to assess the decodability
%   of odor identity from sniffing data using pattern correlation analysis.
%   It computes correlation matrices to compare sniff traces from trials with the same odor versus different odors.
%
%   It also performs SVM classification to decode odor identity from sniffing parameters.
%
% BASIC INPUTS:
%   - nodor: Number of odors (default: 160).
%   - wind: Number of samples (default: 7500). Only when decoding is on raw
%       traces.
%   - dirs: Cell array containing paths to directories with sniff feature data.
%   - behav: Behavioral data loaded from 'NEMO_perceptual2.mat'.
%
% BASIC OUTPUTS:
%   - Figures illustrating:
%       * Pattern correlations between sniffs of the same odor versus different odors.
%       * SVM classification performance.
%   - Variables saved to 'svmpred.mat' containing SVM predictions and performance metrics.
%
% VivekSagar2016@u.northwestern.edu; Nov 27 2024

% To run a demo: supply filepath for mainroot

%% Basic settings
mainroot = 'C:\Work\SFP\Scripts';
rootf = fullfile(mainroot,'supporting_files');
savepath = fullfile(mainroot,'\examples\example_SVM');
nodor = 160;
wind = 7500; % Number of samples
num_subjects = 3;
behav = load(fullfile(rootf,'NEMO_perceptual2.mat'));
corrmat_ = true; % Run featureless patten decoding. Currently not supported at demo version.
sig_clean = true;
%% Basic decoding
if corrmat_
    maindir = 'C:\Work\SFP';
    sig_decorrelate = false;
    corrmoda = zeros(3,1,2); % Subject-wise pattern correlations
    for ss = 1:num_subjects
        subdir = fullfile(maindir,sprintf('sfp_behav_s%02d_correct',ss));
        load(fullfile(subdir,'sfp_feats_main.mat'))
        Fless_mat = vertcat(fless_mat{:});
        % Fless_mat = vertcat(feat_mat{:});

        if ss==3; s2 = 4; else; s2 = ss; end
        onsets = load(fullfile(subdir,sprintf('conditions_NEMO%02d.mat',s2)),'onsets');
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

        if sig_decorrelate
            Fless_mat_pruned = sfp_decorrelate_sig(Fless_mat_pruned,group_vec,behav.behav(ss).ratings(:,1));
            Fless_mat_pruned = sfp_decorrelate_sig(Fless_mat_pruned,group_vec,behav.behav(ss).ratings(:,2));
        end

        % Fless_mat_pruned = Fless_mat(:,[3 4 9:32]);
        Fless_mat_pruned(isnan(Fless_mat_pruned))=0;

        if sig_clean
            Fless_mat_pruned = Fless_mat_pruned-mean(Fless_mat_pruned);
        end

        Fless_corr = corrcoef(Fless_mat_pruned');

        % Trials belong to same odor
        M_on = logical(unity);
        M_on(logical(eye(size(M_on)))) = false; % Trivial correlation

        % Trials belong to different odors
        M_off = ~logical(unity);

        % Correlation across all trials belonging to different odors
        M_mid = Fless_corr;
        M_mid(logical(unity)) = nan;
        M_mid = nanmean(M_mid,2);

        % Correlation across all trials belonging to same odor
        M_mid_on = Fless_corr;
        M_mid_on(~M_on) = nan;
        M_mid_on = nanmean(M_mid_on,2);

        pval = ranksum( M_mid_on,M_mid)  ;
        corrmoda(ss,1,1) = mean(Fless_corr(M_on));
        corrmoda(ss,1,2) = mean(Fless_corr(M_off));
        % corrmodb(ss) = std(M_mid_on-M_mid)/sqrt(4560);
        corrmodp(ss)  = pval;
    end

    % Figure plotting
    figure('Position',[0.5 0.5 320 240])
    rsa_P1 = corrmoda;
    S_mat = squeeze(mean(rsa_P1,1));
    S_err = squeeze(std(rsa_P1,1))./sqrt(3)*1.96;
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
    c_s = {'r','g','b'}; % Data dots for subjects
    for jj = 1:3
        plot([1 2],squeeze(rsa_P1(jj,1,:)),c_s{jj},'handle','off')
    end
    xticks([1 2])
    xticklabels({'Same odor','Different odor'})
    ylabel('Pattern correlation')
end

%% SVM Decoding
nfolds = 10;
numpcs = [13 11 11]; % 90% Variance
svm_trainer2 = true; % Train SVM
rm_trials = false;
% Additional analyses on odor distribution
odorlist = readtable(fullfile(rootf,'Odor_list_new.xlsx'));
odorlist = odorlist(:,1:2);
odorlist_out = fullfile(rootf,'Odor_list_dec.xlsx');
odorlist = rmmissing(odorlist);
M = containers.Map(odorlist.CID,odorlist.Odor);

if svm_trainer2
    corrmod = zeros(3,1);
    predictions = cell(3,1);
    baseline = cell(3,1);
    pvalue = zeros(3,1);
    figure()
    corrmoda = zeros(3,1,2);
    hold on

    tb_odor = {}; % table to store odor information of decoding
    for ss = 1:num_subjects
        subdir = fullfile(rootf,sprintf('sfp_behav_s%02d_correct',ss));
        load(fullfile( subdir,'sfp_feats_main.mat'))

        if ss==3; s2 = 4; else; s2 = ss; end
        onsets_mat = load(fullfile( subdir,sprintf('conditions_NEMO%02d.mat',s2)),'onsets','names');
        onsets = onsets_mat.onsets;
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

        % Fless_mat = vertcat(fless_mat{:});
        % Fless_mat_pruned = Fless_mat(:,1:100:wind);

        Fless_mat = vertcat(feat_mat{:});
        Fless_mat_pruned  = Fless_mat(:,[3 4 9:21 23:31]);
        % Fless_mat_pruned = Fless_mat(:,[3 4 9:31]);

        Fless_mat_pruned(isnan(Fless_mat_pruned))=0;
        Fless_mat_pruned = zscore(Fless_mat_pruned,1);

        [coeff,Fless_mat_pruned,~,~,var] = pca(Fless_mat_pruned);
        Fless_mat_pruned = Fless_mat_pruned(:,1:numpcs(ss));

        oid_ = 1:nodor;
        oid = oid_(group_vec)';

        [~,predictions_vec] = Classify_Permute_VS2(Fless_mat_pruned, oid, 5);

        % behavioral_corr
        accuracies = predictions_vec==oid;

        % distribution of accuracies
        oid_groups = findgroups(oid);
        accuracy_dist = splitapply(@mean, accuracies, oid_groups);
        [~,argsort] = sort( accuracy_dist,'descend');
        odorids = cellfun(@(x) str2num(x), onsets_mat.names,'UniformOutput',false);
        odorids = vertcat(odorids{:});
        odorids = odorids(argsort);
        cell_str = cell(nodor,1); for oid_ = 1:nodor; cell_str{oid_} = M(odorids(oid_)); end
        tb_odor{ss} = table(odorids, cell_str, accuracy_dist(argsort), ...
            'VariableNames', {'CID', 'Odor', 'Accuracy'});
        % writetable( tb_odor{ss}, odorlist_out, 'Sheet', sprintf('Sheet%d', ss));

        % if remove trials with high accuracy
        if rm_trials
            oid_rm = find(accuracy_dist >0.2);
            trial_rm = ismember(group_vec,oid_rm);
            [~,predictions_vec] = Classify_Permute_VS2(Fless_mat_pruned(~trial_rm,:), oid(~trial_rm), 5);
            accuracies = predictions_vec== oid(~trial_rm);
        end

        % Accuracies
        corrmod(ss) = sum(accuracies)/length(accuracies);
        actual_behav = behav.behav(ss).ratings(group_vec,:);
        actual_behav = actual_behav(~accuracies ,1:end);
        predicted_behav = behav.behav(ss).ratings(predictions_vec,:);
        predicted_behav = predicted_behav(~accuracies ,1:end);
        predictions{ss} = iter_corr( actual_behav,predicted_behav);

        % Shuffle test on perceptual correlation on wrong trials
        baseline{ss}=  iter_corr_shuff( actual_behav,predicted_behav);
        corrmoda(ss,1,1) = mean( predictions{ss} );
        corrmoda(ss,1,2)= mean( baseline{ss} ,'all');

        T_shuff = mean(baseline{ss});
        t_stat = mean(predictions{ss});
        pValue_box(ss) = 2 * min(mean(T_shuff >= t_stat), mean(T_shuff <= t_stat));

        subplot(1,3,ss)
        hold on
        histogram(accuracy_dist,10)
        if ss==1
            ylabel('Number odors')
        elseif ss==2
            xlabel('Mean accuracy')
        end
    end
end

% Basic SVM
figure('Position',[0 0 320 240])
hold on
bar(mean(corrmod))
errorbar(mean(corrmod),std(corrmod)./sqrt(3)*1.96)
hold on
c_s = {'r','g','b'};
for ss = 1:3; plot([1],corrmod(ss),c_s{ss},'Marker','.','MarkerSize',15); end
yline(1/160)
ylabel('Performance')
savefig(fullfile(savepath,'svm'))
% print(fullfile(savepath,'svm'),'-dpng')

% Boxplots
figure()
VS_boxplot(predictions)
xticks(1:3)
xticklabels({'S1','S2','S3'})
ylabel('Performance (r)')
savefig('behav')
print('behav','-dpng')
p = [];
for ii = 1:3
    p(ii) = signrank(predictions{ss});
end
savefig(fullfile(savepath,'boxp'))
% print('boxp','-dpng')

% P-value
p_value_svm = arrayfun(@(x) ARC_computePValueOneTailed(x, 160, 4320),corrmod);
p_value_main_svm = ARC_computePValueOneTailed(mean(corrmod), 160, 4320);
SFP_clearLargeVariables
save('svmpred')
