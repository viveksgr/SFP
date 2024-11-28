%--------------------------------------------------------------------------
% SCRIPT: sfp_multilogistic.m
%
% DESCRIPTION:
%   This script performs multilogistic regression analysis to assess how sniffing
%   features predict perceptual ratings of odors. It processes sniffing data and
%   behavioral ratings from multiple subjects, fits logistic regression models,
%   and visualizes the results.
%
% BASIC INPUTS:
%   - Sniffing features data ('sfp_feats_main.mat') for each subject.
%   - Behavioral ratings data ('NEMO_perceptual2.mat').
%
% BASIC OUTPUTS:
%   - Figures displaying logistic regression weights and model performance.
%   - Statistical significance (p-values) of the models.
%
% VivekSagar2016@u.northwestern.edu; Nov 27 2024
% To run a demo: supply filepath for mainroot

%% Multilogistic regression

% Sniffing modulation
grp_ = true; % Runs multilogistic model
num_descrip = 31;
numsubjects = 3;
mainroot = 'C:\Work\SFP\Scripts';
rootf = fullfile(mainroot,'supporting_files');
savepath = fullfile(mainroot,'\examples\example_multilogistic2');
mkdir(savepath)

settings_.pcamaker = true;
numpcs = [13 11 11]; % 90% Variance from each subject

figure()
hold on
load(fullfile(rootf,'snifflabels.mat'))
proper_list(16)=[]; % This descriptor is removed
if grp_
    nodor = 160;
    wind = 3500; % Number of samples
    color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
    behav = load(fullfile(rootf,'NEMO_perceptual2.mat'));
    ndisc = size(behav.behav(1).ratings,2);

    num_eff = zeros(3,ndisc);
    num_pval = zeros(3,ndisc);
    for ss = 1:numsubjects
        fprintf('subject %02d: \n',ss)
        subdir = fullfile(rootf,sprintf('sfp_behav_s%02d_correct',ss));
        load(fullfile( subdir,'sfp_feats_main.mat'))
        % Fless_mat = vertcat(fless_mat{:});
        feat_mat = vertcat(feat_mat{:});
        feat_mat_pruned = feat_mat(:,[3 4 9:21 23:num_descrip]);
        feat_mat_pruned(isnan(feat_mat_pruned))=0;
        feat_mat_pruned = zscore(feat_mat_pruned);


        if ss==3; s2 = 4; else; s2 = ss; end
        onsets = load(fullfile( subdir,sprintf('conditions_NEMO%02d.mat',s2)),'onsets');
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

        nfeatures = size(feat_mat_pruned,2);
        if settings_.pcamaker
            npc = numpcs(ss);
           
            [coeff,sc,~,~,var] = pca(feat_mat_pruned);
            feat_mat_pruned = sc(:,1:npc);

            % varsum = cumsum(var);
            % varmat{ss} = varsum;
            % vardiff{ss} = var;
            % subplot(1,3,ss)
            % hold on
            % plot(varsum)
            % xlabel('Num descriptors')
            % ylabel('PCA variance explained')
        end

        wt_sc_mat = zeros(nfeatures,ndisc);
        for pp = 1: ndisc      
            behav_ratings = behav.behav(ss).ratings(:,pp);
            behav_ratings_= behav_ratings(group_vec);

            nbins = 3;
            delta = 0.0001*linspace(0,1,nbins+1);
            edges_ = quantile(behav_ratings_, linspace(0, 1, nbins + 1));
            edges_ = edges_+delta;
            y = discretize(behav_ratings_, edges_);


            [num_eff(ss,pp) , num_pval(ss,pp) ,wt] = SFP_logisticRegressionCV_ordinal(y, feat_mat_pruned);  
            
            if settings_.pcamaker
                coeffmat = coeff(:,1:npc);
                wt_sc = coeffmat.*wt';
                wt_sc_mat(:,pp) = abs(mean(wt_sc,2));
            else
                wt_sc_mat(:,pp) = wt;
            end
        end

        wt_mat(ss,:) = mean(wt_sc_mat,2);
        wt_mat_cell{ss} = wt_sc_mat;

        % subplot(1,3,ss)
        % hold on
        % axis tight
        % imagesc(wt_sc_mat)
        % if ss==1
        % yticks(1:1: nfeatures)
        % yticklabels(strrep(proper_list, '_', ' '))
        % end
        % xticks(1:1:ndisc   )
        % xtickangle(90)
        % xticklabels(behav.behav(ss).percepts)

    end
end 

nS = 3;
behav_lab1 = behav.behav(1).percepts;
behav_lab2 = behav.behav(2).percepts;

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
rels_mean = nanmean(rels,1);
[~,argsort] = sort(rels_mean,'descend');
behav_labs_sort = behav_labs(argsort);
rels_sort = rels(:,argsort);

figure()
bar(1:21,nanmean(rels_sort))
hold on
errorbar(1:21,nanmean(rels_sort),1.96*nanstd(rels_sort)/sqrt(3),'.')
c = {'r.','g.','b.'};
for ss = 1:3
    plot(1:21,rels_sort(ss,:),c{ss})
end
xticks(1:21)
xticklabels(behav_labs_sort)
xtickangle(90)
yline(1/3)
savefig(fullfile(savepath,'logistic.fig'))

% FDR correction
num_pval_cell = mat2cell(num_pval,[1 1 1]);
[a,num_p_fdr] = cellfun(@fdr_benjhoc,num_pval_cell,'UniformOutput',false);
