%% RSA time-series
corrmat_ = true;
dwnsample = 100;
if corrmat_
    nodor = 160;
    wind = 75; % Number of samples
    dirs = {'C:\Work\SFP\sfp_behav_s01_correct';
        'C:\Work\SFP\sfp_behav_s02_correct';
        'C:\Work\SFP\sfp_behav_s04_correct'};
    %   color = [0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330; 0 0.4470 0.7410];
    behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));

    hold on
    corrmod = zeros(3,wind);
    corrmod_int = zeros(3,wind);
    corrmod_pls = zeros(3,wind);
    sniff_trace = zeros(3,wind);
    for ss = 1:length(dirs)
        load(fullfile(dirs{ss},'sfp_feats_main.mat'))
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

        % Behavioral features
        behav_ratings = behav.behav(ss).ratings;
        behav_ratings_ =  behav_ratings(group_vec,:);
        b1 = corrcoef(behav_ratings_')';

        rint = regressmeout(behav_ratings(:,[2:end])',repmat(behav_ratings(:,1)',17,1));
        rint = rint(:,group_vec);
        brint = corrcoef(rint);

        rpls = regressmeout(behav_ratings(:,[1 3:end])',repmat(behav_ratings(:,2)',17,1));
        rpls = rpls(:,group_vec);
        bpls = corrcoef(rpls);

        for jj = 1:wind
            Fless_corr = -abs(Fless_mat_pruned(:,jj)-Fless_mat_pruned(:,jj)');
            b1_vec = b1(utl_mask);
            b1_vec = regressmeout(b1_vec',unity(utl_mask)')';

            temp = corrcoef(b1_vec ,Fless_corr( utl_mask ));
            corrmod(ss,jj) = temp(2);

            % temp = corrcoef(brint,Fless_corr);
            % corrmod_int(ss,jj) = temp(2);

            % temp = corrcoef(bpls,Fless_corr);
            % corrmod_pls(ss,jj) = temp(2);
        end
        sniff_trace(ss,:) = mean(Fless_mat_pruned ,1);
    end
end

figure('Position',[0.5 0.5 640 480])
hold on
yyaxis left
plot((1:wind)/10,mean(sniff_trace,1))
ylabel('Sniff Trace')

yyaxis right
plot((1:wind)/10,mean(corrmod,1))
% plot((1:wind)/10,mean(corrmod_int,1))
% plot((1:wind)/10,mean(corrmod_pls,1))
ylabel('RSA performance')

legend({'Sniff','Perceptual Similarity'})
yyaxis left
xlabel('time(s)')

figure()
hold on
for ii = 1:3
    subplot(1,3,ii)
    hold on
    yyaxis left
    plot((1:wind)/10,corrmod(ii,:))
    plot((1:wind)/10,corrmod_int(ii,:))
    plot((1:wind)/10,corrmod_pls(ii,:))
    ylabel('RSA performance')
    yyaxis right
    plot((1:wind)/10,sniff_trace(ii,:))
    ylabel('Sniff Trace')
    legend({'Basic','-Int','-Pls','Sniff'})
    yyaxis left
    title(sprintf('Subject: %02d',ii))
end

figure()
hold on
x = mean(sniff_trace,1);
y = mean(corrmod,1);
c = linspace(1,10,wind);
scatter(x(1:75),y(1:75),25,c(1:75))
ylabel('Rep Similarity')
xlabel('Sniff')
% 
figure('Position',[0 0 720 200])
hold on
for ii = 1:3
    subplot(1,3,ii)
    hold on
    x = sniff_trace(ii,:);
    y = corrmod(ii,:);
    c = linspace(1,10,wind);
    scatter(x,y,25,c)
    title(sprintf('Subject: %02d',ii))
end
savefig('pls_sub_cross')
print('pls_sub_cross','-dpng')