%% RSA time-series
corrmat_ = true;
dwnsample = 100;
switcher = 'basic'; %'basic', 'int' or 'pls';

if corrmat_
    nodor = 160;
    wind = 75; % Number of samples
    t_idx = 1:dwnsample:7500;
    win = 50;
    dirs = {'C:\Work\SFP\sfp_behav_s01_correct';
        'C:\Work\SFP\sfp_behav_s02_correct';
        'C:\Work\SFP\sfp_behav_s04_correct'};
    %   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
    behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));

    hold on
    corrmod = zeros(3,wind);
    corrmod_t = zeros(3,wind);
    corrmodu = zeros(3,wind);
    corrmodu_t = zeros(3,wind);
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
          % behav_ratings_  = zscore(  behav_ratings_ ,[],1);

        switch switcher
            case 'basic'
                b1 = corrcoef(behav_ratings_')';
                % b1 = pdist(behav_ratings_);
                % b1 = squareform(b1);
                DM_mat = [b1(utl_mask) unity(utl_mask) task_run(utl_mask) sess_run(utl_mask) set_run(utl_mask)];
            case 'int'
                rint = regressmeout(behav_ratings(:,[2:end])',repmat(behav_ratings(:,1)',17,1));
                rint = rint(:,group_vec);
                b1 = corrcoef(rint);
            case 'pls'
                rpls = regressmeout(behav_ratings(:,[1 3:end])',repmat(behav_ratings(:,2)',17,1));
                rpls = rpls(:,group_vec);
                b1 = corrcoef(rpls);
            otherwise
                error('Wrong condition')
        end
        fprintf('Subject: %02d',ss)
        fprintf('\n')
        for jj=1:wind; fprintf('.'); end
        fprintf('\n')
        for jj = 1:wind
            fprintf('|')
            % Fless_corr = -abs(Fless_mat_pruned(:,jj)-Fless_mat_pruned(:,jj)');
             % Fless_corr = corrcoef(Fless_mat(:,t_idx(jj):t_idx(jj)+win)');

            Fless_corr = -pdist(Fless_mat(:,t_idx(jj):t_idx(jj)+win));
            Fless_corr = squareform(Fless_corr);

            [wt2,t_sc2] = ARC_multicomputeWeights_tsc(DM_mat, Fless_corr(utl_mask));
            corrmod(ss,jj) = wt2(2);
            corrmod_t(ss,jj) = t_sc2(2);

            corrmodu(ss,jj) = wt2(3);
            corrmodu_t(ss,jj) = t_sc2(3);
        end
        sniff_trace(ss,:) = mean(Fless_mat_pruned ,1);
        fprintf('\n')
    end
end

figure('Position',[0.5 0.5 1280 320])
hold on
nsig = [nchoosek(4560,2),nchoosek(4320,2),nchoosek(4320,2)];
for ii = 1:3
    subplot(1,3,ii)
    hold on
    yyaxis right
    plot((1:wind)/10,corrmod(ii,:),'LineWidth',0.2)
    t_thr = tinv(0.975,nsig(ii));
    sig_t = corrmod_t(ii,:)>t_thr;
    main_t = corrmod(ii,:);
    main_t(~sig_t)=nan;
    plot((1:wind)/10,main_t,'Linestyle',"-",'LineWidth',2)

    % plot((1:wind)/10,corrmodu(ii,:),'g','LineWidth',0.2)
    % t_thr = tinv(0.975,nsig(ii));
    % sig_t = corrmodu_t(ii,:)>t_thr;
    % main_t = corrmodu(ii,:);
    % main_t(~sig_t)=nan;
    % plot((1:wind)/10,main_t,'Linestyle','g',"-",'LineWidth',2)

    ylabel('RSA performance')
    yyaxis left
    plot((1:wind)/10,sniff_trace(ii,:))
    ylabel('Sniff Trace')
    yyaxis left
    title(sprintf('Subject: %02d',ii))
end
savefig('RSA_maintrace_subwise')
print('RSA_maintrace_subwise','-dpng')
print(gcf,'-vector','-dsvg',[fullfile(pwd,'RSA_maintrace_subwise'),'.svg']) % svg

figure('Position',[0.5 0.5 640 480])
[corrmod_t_adj, ~,~] = fdr_t_scores(corrmod_t,[sum(utl_mask(:)) sum(utl_mask(:)) sum(utl_mask(:))] );
SFP_resetfigs
hold on
yyaxis left
plot((1:wind)/10,mean(sniff_trace,1))
ylabel('Sniff Trace')
yyaxis right

plot((1:wind)/10,mean(corrmod,1))
sig_t = corrmod_t_adj>t_thr;
sig_t = sum(sig_t,1)>=2;
main_t = mean(corrmod,1);
main_t(~sig_t)=nan;
plot((1:wind)/10,main_t,'Linestyle',"-",'LineWidth',2)

% plot((1:wind)/10,mean(corrmodu,1),'g')
% sig_t = corrmodu_t>t_thr;
% sig_t = sum(sig_t,1)>=2;
% main_t = mean(corrmodu,1);
% main_t(~sig_t)=nan;
% plot((1:wind)/10,main_t,'g','Linestyle',"-",'LineWidth',2)

ylabel('RSA performance')
legend({'Sniff','Perceptual Similarity'})
yyaxis left
xlabel('time(s)')
savefig('RSA_maintrace')
print('RSA_maintrace','-dpng')
print(gcf,'-vector','-dsvg',[fullfile(pwd,'RSA_maintrace'),'.svg']) % svg

% % Plots2
% figure('Position', [0.5 0.5 640 480]);
% ax1 = axes;
% hold(ax1, 'on');
% plot(ax1, (1:wind)/10, mean(sniff_trace, 1), 'b');
% ylabel(ax1, 'Sniff Trace');
% ax2 = axes('Position', ax1.Position, 'Color', 'none', 'YAxisLocation', 'right', 'XTick', []);
% hold(ax2, 'on');
% plot(ax2, (1:wind)/10, mean(corrmod, 1), 'r');
% ylabel(ax2, 'Correlation Model');
% xlabel(ax1, 'time(s)');
% linkaxes([ax1, ax2], 'x');
% set([ax1, ax2], 'FontName', 'Arial');
% print(gcf, '-vector', '-dsvg', fullfile(pwd, 'RSA_maintrace2.svg'));
% 
% 
% figure()
% hold on
% x = mean(sniff_trace,1);
% y = mean(corrmod,1);
% c = linspace(1,10,wind);
% % plot(x(1:75),y(1:75))
% scatter(x(1:75),y(1:75),25,c(1:75))
% ylabel('Rep Similarity')
% xlabel('Sniff')
% savefig('cross')
% print('cross','-dpng')
% % plot2svg(fullfile(pwd,'cross.svg'))
% print(gcf,'-vector','-dsvg',[fullfile(pwd,'cross'),'.svg']) % svg
% % 
% figure('Position',[0 0 720 200])
% hold on
% for ii = 1:3
%     subplot(1,3,ii)
%     hold on
%     x = sniff_trace(ii,:);
%     y = corrmod(ii,:);
%     c = linspace(1,10,wind);
%     scatter(x,y,25,c)
%     title(sprintf('Subject: %02d',ii))
% end
% savefig('cross_subwise')
% print('cross_subwise','-dpng')
% print(gcf,'-vector','-dsvg',[fullfile(pwd,'cross_subwise'),'.svg']) % svg
SFP_clearLargeVariables
save('mainmat')