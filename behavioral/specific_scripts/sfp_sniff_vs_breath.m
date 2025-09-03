cd('F:\NEMO_oldrepo\Older_repo-II\NEMO\NEMO_02\breathing\imaging_task\')
kk = 0;
corr_mat = [];
corr_val = [];
num_tr = [];
num_samp = [];
for set = [2:4];
for ii = 1:4
    kk = kk+1;
    datname =fullfile(sprintf('set_%02d',set),'sess_01',sprintf('NEMO02_set%02d_sess01_run%02d.mat',set,ii));
    RI_data = ReadLabChartMat(datname);

    R = RI_data.data{4};
    R = smoothdata(R,'movmean',100)';
    K.RT = 1/1000;
    K.row = ones(length(R),1);
    K.HParam = 50;
    R = spm_filter(K,R);
    R = zscore(R);

    % Breathing belt
    R_diff = diff(R);
    R_t= [R_diff(1); R_diff];
    bmObj = breathmetrics(R_t, 1000, 'humanAirflow');
    bmObj.estimateAllFeatures();
    prop_list = properties(bmObj);
    proper_list = prop_list(7:24);
    feat_superset = [];
    for feat = 1:length(proper_list)
        eval(sprintf('feat_superset = cat(1,feat_superset,bmObj.%s);', ...
            proper_list{feat}))
    end
    feat_superset_br =  feat_superset([3 4 9:14],:)';
     % Spirometer belt
    onsets = feat_superset(1,:)';

        S = RI_data.data{6};
          % Breathing belt
 
    bmObj = breathmetrics(S, 1000, 'humanAirflow');
    bmObj.estimateAllFeatures();
    prop_list = properties(bmObj);
    proper_list = prop_list(7:24);
    feat_superset = [];
    for feat = 1:length(proper_list)
        eval(sprintf('feat_superset = cat(1,feat_superset,bmObj.%s);', ...
            proper_list{feat}))
    end
    t_star = find_nearest(feat_superset(1,:)',onsets); 


    feat_superset_spr =  feat_superset([3 4 9:14], t_star)';
    corr_mat(:,kk) = iter_corr( feat_superset_br',feat_superset_spr');
    corr_val(kk) = fastcorr( S,R_t);
    num_tr(kk) = length(feat_superset_spr);
    num_samp(kk) = length(S);

end
end

figure()
hold on
num_trials = length(corr_val);
bar(mean(corr_val))
errorbar(mean(corr_val),1.96*(std(corr_val)./sqrt(num_trials)))
hold on
for ss = 1:num_trials; plot([1],corr_val(ss),'k','Marker','.','MarkerSize',15); end 
yline(r2t(0.025,min(num_samp)))
ylabel('Performance')
savefig(fullfile(savepath,'svm'))

% % figure()
% R = smoothdata(R,'movmean',2500)';
% K.RT = 1/1000;
% K.row = ones(length(R),1);
% K.HParam = 50;
% R = spm_filter(K,R);
% R = zscore(R);
% plot(R)
% 
% plot(S)
figure()
hold on
yyaxis left
plot([0; diff(R)])
yyaxis right
plot(S)

figure()
bar(mean(corr_mat'))
hold on
errorbar(mean(corr_mat'),std(corr_mat')/sqrt(12),'.')
for zz = 1:8
    plot(ones(1,12)*zz,corr_mat(zz,:),'.k')
end
yline(r2t(0.001,mean(num_tr)))
xticks(1:8)
xtickangle(45)
xticklabels(proper_list([3:4 9:14]))

%% Compare features
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

