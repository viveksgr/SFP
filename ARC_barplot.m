function ARC_barplot(rsa_P1,sw)


if nargin<2
    sw = true;
end

S_mat = squeeze(nanmean(rsa_P1));
S_err = squeeze(nanstd(rsa_P1))./sqrt(3);
if sw
figure('Position',[0.5 0.5 400 250])
hold on
end
ngroups = size(S_mat, 1);
nbars = size(S_mat, 2);
bar(S_mat);
colororder("earth")
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
x_m = [];
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/(2*nbars);
    errorbar(x, S_mat(:,i), S_err(:,i), 'k.');
    x_m = [x_m; x];
end

% legend({'Perceptual','Chemical','Mutual'})
% legend()
% Subject data points
c_s = {'r','g','b'}; % Data dots for subjects
for ii = 1:size(rsa_P1,2) % For bars for perceptual, chemical and combinations
    for jj = 1:3
        plot(x_m(:,ii),squeeze(rsa_P1(jj,ii,:)),c_s{jj},'handle','off')
    end
end

% yline(r2t(0.05,sum(utl_mask2(:))));